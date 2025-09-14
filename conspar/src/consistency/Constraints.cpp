/**
 * Author: Mario Joao Jr.
 * 
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <iostream>
#include <fstream>
#include <exception>
#include <cmath>
#include <sstream>
#include "Constraints.h"

using namespace std;

namespace clustalw
{

Constraints::Constraints(int nSeqs)
{
	vector<vector<SparseMatrix *> > sparseMat (nSeqs, vector<SparseMatrix *>(nSeqs, NULL));
	sparseMatrices = sparseMat;
	this->nSeqs = nSeqs;
	this->nEdges = 0;
	// resIndex is a matrix indicating where in constrList s1,s2 list of constraints starts.
	this->resLocks.resize (nSeqs);
	this->resIndex.resize (nSeqs);
	this->resEdges.resize (nSeqs);
	for (int i=0; i < nSeqs; i++)
	{
		this->resIndex[i].resize (nSeqs, -1);
		this->resEdges[i] = {};
		omp_init_lock (&(this->resLocks[i]));
	}
        this->maxValue = 0;
        this->maxExtValue = 0;
        this->noMatch = 0;
        this->normalise = 1000;
        this->extEdges = new vector<ExtEdge>;
}

int Constraints::getIniCLPair (int s1, int s2)
{
	if (resIndex[s1][s2] != -1)
		return resIndex[s1][s2];

	for (int i = s1; i < nSeqs; i++)
	{
		for (int j = s2; j < nSeqs; j++)
		{
			if (resIndex[i][j] != -1)
				return resIndex[i][j];
		}
	}
		
	return constrList.size();
}

void Constraints::moveRight (int ini)
{
	ConstrEntry     *entry;

	entry = &constrList[constrList.size() - 1];
//	entry->setInd(entry->getInd()+1);
	constrList.push_back (*entry);

	for (int i=constrList.size() - 1; i > ini; i--)
	{
		constrList[i] = constrList[i-1];
		constrList[i].setInd(i);
	}
} 

void Constraints::addEntry2 (ConstrEntry entry)
{
	//adds an entry to the list
	int	iCLPair,
		i;

	if (entry.getInd() != -1)
	{
		constrList[entry.getInd()] = entry;
		return;
	}

	iCLPair = getIniCLPair(entry.getS1(), entry.getS2());

//	printf ("iCLPair = %d, constrList.size=%lu\n", iCLPair, constrList.size());

	if (iCLPair == constrList.size())
	{
		resIndex[entry.getS1()][entry.getS2()] = iCLPair;
		entry.setInd(iCLPair);
		constrList.push_back(entry);
		return;
	}

	for (i = iCLPair; (i < constrList.size()) && (constrList[i].getS1() < entry.getS1()); i++);

	if (constrList[i].getS1() == entry.getS1())
	{
		for (; (i < constrList.size()) && (constrList[i].getS2() < entry.getS2()); i++);

		if (constrList[i].getS2() == entry.getS2())
		{
			for (; (i < constrList.size()) && (constrList[i].getR1() < entry.getR1()); i++);

			if (constrList[i].getR1() == entry.getR1())
			{
				for (; (i < constrList.size()) && (constrList[i].getR2() <= entry.getR2()); i++);
			}
		}
	}

/*
	for (i = iCLPair; (i < constrList.size()); i++)
	{
		if (constrList[i].getS1() > entry.getS1())
			break;
		if (constrList[i].getS2() > entry.getS2())
			break;
		if (constrList[i].getR1() > entry.getR1())
			break;
		if (constrList[i].getR2() > entry.getR2())
			break;
	}

	printf ("Moveright(%d)\n", i);
*/

	moveRight (i);
	entry.setInd(i);
	constrList[i] = entry;

/*
  if (i==1){CL=insert_entry (entry, CL,1);}
  else
    {
      i=get_entry_index(entry,CL);
      if (i<0)insert_entry (entry,CL,-i);
      else update_entry (entry, CL, i);
    }
*/
}

void Constraints::addEntry (ConstrEntry entry)
{
//	ConstrEntry	entry2 (entry.getS2(), entry.getS1(), entry.getR2(), entry.getR1(), entry.getWeight(), entry.getCons(), -1);
	//adds an entry and its mirror to the list
	//if INDEX is set the entry replaces the entry with a similar index
	//otherwise the entry (and its mirror) are added
//  if (entry[INDEX])return add_entry2list2(entry, CL);

	if (entry.getInd() != -1)
	{
		this->constrList[entry.getInd()] = entry;
		this->setEdge(resEdges[entry.getS1()][entry.getR1()], entry.getS2(), entry.getR2(), entry.getWeight());
		this->setEdge(resEdges[entry.getS2()][entry.getR2()], entry.getS1(), entry.getR1(), entry.getWeight());
		return;
	}

	// Colocar duas entradas ou apenas colocar os índices apontando para o mesmo lugar?

//	this->addEntry2 (entry2);

	/* 
	 * A constrList é inutil. Toda a estrutura está representada por Edges. MJ 04/11/2023.
	 */
//	this->addEntry2 (entry);
	this->addEdge (entry);
}

void Constraints::setEdge (vector<Edge> &vet, int s2, int r2, int weight)
{
	for (int i = 0; i < vet.size(); i++)
	{
		if ((vet[i].getS2() == s2) && (vet[i].getR2() == r2)) 
		{
			vet[i].setEdge (s2, r2, weight);
		}
	}
}

void Constraints::addEdge (ConstrEntry entry)
{
	int	s1 = entry.getS1(),
	 	s2 = entry.getS2(),
	 	r1 = entry.getR1(),
	 	r2 = entry.getR2(),
	 	weight = entry.getWeight(),
		i;
	int	len;

//	omp_set_lock (&(resLocks[s1]));
//	omp_set_lock (&(resLocks[0]));
	#pragma omp critical
	{
	len = resEdges[s1].size();
	if (r1 + 1 > len)
	{
		resEdges[s1].resize (r1 + 1);
	}
	Edge e1 (s2, r2, weight);
//	resEdges[s1][r1].push_back(e1);

	resEdges[s1][r1].resize (resEdges[s1][r1].size () + 1);
	for (i = resEdges[s1][r1].size() - 1; i > 0; i--)
	{
		if (resEdges[s1][r1][i-1].getS2() <  e1.getS2())
			break;

		if ((resEdges[s1][r1][i-1].getS2() ==  e1.getS2()) && (resEdges[s1][r1][i-1].getR2() <  e1.getR2()))
			break;

		resEdges[s1][r1][i] = resEdges[s1][r1][i-1];
	}

	resEdges[s1][r1][i] = e1;
//	omp_unset_lock (&(resLocks[s1]));

//	omp_set_lock (&(resLocks[s2]));
	len = resEdges[s2].size();
	if (r2 + 1 > resEdges[s2].size())
	{
		resEdges[s2].resize (r2 + 1);
	}
	Edge e2 (s1, r1, weight);
//	resEdges[s2][r2].push_back(e2);

	resEdges[s2][r2].resize (resEdges[s2][r2].size () + 1);
	for (i = resEdges[s2][r2].size() - 1; i > 0; i--)
	{
		if (resEdges[s2][r2][i-1].getS2() <  e2.getS2())
			break;

		if ((resEdges[s2][r2][i-1].getS2() ==  e2.getS2()) && (resEdges[s2][r2][i-1].getR2() <  e2.getR2()))
			break;

		resEdges[s2][r2][i] = resEdges[s2][r2][i-1];
	}

	resEdges[s2][r2][i] = e2;
	}
//	omp_unset_lock (&(resLocks[0]));
//	omp_unset_lock (&(resLocks[s2]));

/*
//	#pragma omp critical
	{
	nEdges += 2;
	}
*/
}

void Constraints::saveEdges (string fileName)
{
	ofstream outFile;

	try
	{
		outFile.open (fileName, ios::out | ios::trunc);
		outFile.setf(ios::fixed);

		for (int i = 0; i < resEdges.size(); i++)
		{
			for (int j = 0; j < resEdges[i].size(); j++)
			{
				if (!resEdges[i][j].size())
					continue;

				for (int k = 0; k < resEdges[i][j].size(); k++)
				{
					Edge	e = resEdges[i][j][k];
					outFile << i << " " << j << " " << e.getS2() << " " << e.getR2() << " " << e.getWeight() << "\n";
				}
			}
		}
	}
	catch (const ofstream::failure& e)
	{
		cout << "Exception writing to Contraint List file";
	}

}

void Constraints::calcNEdges (void)
{
	nEdges = 0;
	for (int i = 0; i < resEdges.size(); i++)
	{
		for (int j = 0; j < resEdges[i].size(); j++)
		{
			nEdges += resEdges[i][j].size();
		}
	}
}

/*
 *	For debug purposes only.
 */
void Constraints::printEdges (void)
{
	int	count = 0;
	printf ("\n");
	for (int i = 0; i < resEdges.size(); i++)
	{
		for (int j = 0; j < resEdges[i].size(); j++)
		{
			if (!resEdges[i][j].size())
				continue;

			printf ("S=%03d,R=%03d==>", i, j);
			for (int k = 0; k < resEdges[i][j].size(); k++)
			{
				Edge	e = resEdges[i][j][k];
				printf ("(%d %d %d) ", e.getS2(), e.getR2(), e.getWeight());
				count++;
			}
			printf ("\n");
		}
	}
	printf ("Number of Edges: %d\n", count);
}

void Constraints::print (void)
{
	printf ("\n");
	for (int i = 0; i < constrList.size(); i++)
	{
		printf ("(%d %d %d %d %d %d %d) ", constrList[i].getS1(), constrList[i].getS2(), constrList[i].getR1(), constrList[i].getR2(), constrList[i].getWeight(), constrList[i].getCons(), constrList[i].getInd());
	}
	printf ("\n");

}

/*
 *	For debug purposes only.
 */
void Constraints::printResIndex (void)
{
	for (int i = 0; i < nSeqs; i++)
	{
		for (int j = 0; j < nSeqs; j++)
		{
			printf ("%d\t", resIndex[i][j]);
		}
		printf ("\n");
	}
}

}

