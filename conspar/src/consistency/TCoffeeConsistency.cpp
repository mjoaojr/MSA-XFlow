/**
 * Author: Mario Joao Jr.
 * 
 * Adapted from TCoffee Code
 */

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <stdlib.h>
#include <omp.h>

#include "TCoffeeConsistency.h"
#include "../alignment/Alignment.h"
#include "Constraints.h"


#define ind(x,y)	(((x)*(maxLen+1))+(y))


void relaxCL (clustalw::Alignment *alignPtr)
{
	int					_numSeqs = alignPtr->getNumSeqs(),
						nEdges,
						finalEdges,
						maxLen;
	clustalw::Constraints			*pCL;
	vector<vector<vector<clustalw::Edge>>>	*pEdges;
//	vector<vector<int>>			hasch;
	vector<vector<int> *>			scores;

	pCL = alignPtr->getPCL();
	pEdges = pCL->getEdges();
/*
	nEdges = pCL->getNEdges();
	printf ("nEdges = %d\n", nEdges);
*/
	pCL->calcNEdges();
	nEdges = pCL->getNEdges();
//	printf ("nEdges = %d\n", nEdges);
	maxLen = alignPtr->getLengthLongestSequence();

/*
	hasch.resize (_numSeqs);
	for (int s1 = 0; s1 < _numSeqs; s1++)
	{
		hasch[s1].resize (maxLen+1);
	}
*/
	scores.resize (_numSeqs);
	for (int s1 = 0; s1 < _numSeqs; s1++)
	{
		scores[s1] = new vector<int>;
	}
	

	omp_set_num_threads(clustalw::userParameters->getNthr());
//	finalEdges = 0;
	#pragma omp parallel for schedule (dynamic, 1)
//	#pragma omp parallel for schedule (static, 1)
//	#pragma omp parallel for schedule (static)
//	#pragma omp parallel for schedule (guided, 1)
	for (int s1 = 0; s1 < _numSeqs; s1++)
	{
		int	*hasch;
		vector<int>	*pscores;

		hasch = (int *) calloc (_numSeqs*(maxLen+1), sizeof(int));
		pscores = scores[s1];
		for (int r1 = 1; r1 <= alignPtr->getSeqLength (s1 + 1) && r1 < (*pEdges)[s1].size(); r1++)
		{
			int	len1,
				norm1;
			norm1=0;
			len1 = (*pEdges)[s1][r1].size();
			for (int x = 0; x < len1; x++)
			{
				int	t_s1,
					t_r1;
				// Não faz sentido hasch[0][27] ser = 4
				t_s1 = (*pEdges)[s1][r1][x].getS2();
				t_r1 = (*pEdges)[s1][r1][x].getR2();
//				hasch[t_s1][t_r1] = (*pEdges)[s1][r1][x].getWeight();
				hasch[ind(t_s1,t_r1)] = (*pEdges)[s1][r1][x].getWeight();
				norm1++;
			}

			for (int a = 0; a < len1; a++)
			{
				int	s2,
					r2,
					len2,
					norm,
					norm2;
				double	score;
				score = 0.0;
				norm2 = 0;
				s2 = (*pEdges)[s1][r1][a].getS2();
				r2 = (*pEdges)[s1][r1][a].getR2();
				len2 = (*pEdges)[s2][r2].size();

				for (int x = 0; x < len2; x++)
				{
					int	t_s2,
						t_r2;
					t_s2 = (*pEdges)[s2][r2][x].getS2();
					t_r2 = (*pEdges)[s2][r2][x].getR2();

					if (t_s2 == s1 && t_r2 == r1)
					{
						score += (float)(*pEdges)[s2][r2][x].getWeight();
					}
					else 
					{
//						if (hasch[t_s2][t_r2])
						if (hasch[ind(t_s2,t_r2)])
						{
							score += min(((float)hasch[ind(t_s2,t_r2)]),((float)(*pEdges)[s2][r2][x].getWeight()));
						}
					}
					norm2++;
				}

				norm = min (norm1, norm2);
				/*
				if (normalisation_mode==0)
				norm=MIN(norm1,norm2);
				else if (normalisation_mode==1)
				norm=MAX(norm1,norm2);
				else if (normalisation_mode==2)
				norm=(CL->S)->nseq;
				*/

				score = ((norm) ? (score / norm) : 0);
				pscores->push_back ((int)score);
//				scores[s1].push_back ((int)score);
//				finalEdges++;
//				printf ("%d ", (int)(score));
			}

			for (int x = 0; x < len1; x++)
			{
				int	t_s1,
					t_r1;
				t_s1 = (*pEdges)[s1][r1][x].getS2();
				t_r1 = (*pEdges)[s1][r1][x].getR2();
//				hasch[t_s1][t_r1] = 0;
				hasch[ind(t_s1,t_r1)] = 0;
			}
		}

		free (hasch);
	}

//	printf ("\nnEdges = %d e finalEdges = %d\n", nEdges, finalEdges);


	/*
	 * Updates Edges and filter it at same time. Threshold = 10. It is hardcoded.
	 * Analyse if it is good to use it hardcoded or used a parameter.
	 * Does this threshold influences the result?
	 */
	finalEdges = nEdges;
	#pragma omp parallel for schedule (dynamic, 1)
//	#pragma omp parallel for schedule (static, 1)
//	#pragma omp parallel for schedule (static)
//	#pragma omp parallel for schedule (guided, 1)
	for (int s1 = 0; s1 < _numSeqs; s1++)
	{
		int	a = 0,
			difedges = 0;
		vector<int>	*pscores;

		pscores = scores[s1];
		for (int r1 = 0; r1 <= alignPtr->getSeqLength (s1 + 1) && r1 < (*pEdges)[s1].size(); r1++)
		{
			int	len1;
			len1 = (*pEdges)[s1][r1].size();
			for (int x = 0; x < len1;)
			{
//				t_s1 = (*pEdges)[s1][r1][x].getS2();
//				t_r1 = (*pEdges)[s1][r1][x].getR2();
//				printf ("%d %d %d %d %d\n", s1, r1, t_s1, t_r1, scores[a]);
//				if (scores[s1][a] > 10)
				if ((*pscores)[a] > 10)
				{
//					(*pEdges)[s1][r1][x].setWeight(scores[s1][a++]);
					(*pEdges)[s1][r1][x].setWeight((*pscores)[a++]);
					x++;
				}
				else
				{
					(*pEdges)[s1][r1].erase ((*pEdges)[s1][r1].begin() + x);
//					finalEdges--;
					difedges++;
					a++;
					len1--;
				}
			}
		}

		#pragma omp critical
		{
		finalEdges -= difedges;
		}
	}

	printf ("\nRelaxation Summary: [%d]--->[%d]\n", nEdges, finalEdges);
	pCL->setNEdges (finalEdges);


	/*
	 * Free hasch memory
	for (int s1 = 0; s1 < _numSeqs; s1++)
	{
		vector<int>().swap(hasch[s1]);
	}
	vector<vector<int>>().swap(hasch);
	 */
	for (int s1 = 0; s1 < _numSeqs; s1++)
	{
		vector<int>().swap(*(scores[s1]));
	}
	vector<vector<int>*>().swap(scores);
	
}

void maxConstraints (clustalw::Alignment *alignPtr)
{
	int					s1,
						r1,
						s2,
						r2,
						w2,
						len1,
						x,
						_numSeqs = alignPtr->getNumSeqs(),
						nEdges,
						maxLen;
	double					score;
	clustalw::Constraints			*pCL;
	vector<vector<vector<clustalw::Edge>>>	*pEdges;
	vector<vector<int>>			maxRes;
	vector<int>				scores;

	pCL = alignPtr->getPCL();
	pEdges = pCL->getEdges();
//	nEdges = pCL->getNEdges();
	maxLen = alignPtr->getLengthLongestSequence();

	maxRes.resize (_numSeqs);
	for (s1 = 0; s1 < _numSeqs; s1++)
	{
		maxRes[s1].resize (maxLen+1);
	}
	

	for (s1 = 0; s1 < _numSeqs; s1++)
	{
		for (r1 = 1; r1 <= alignPtr->getSeqLength (s1 + 1) && r1 < (*pEdges)[s1].size(); r1++)
		{
			len1 = (*pEdges)[s1][r1].size();
			for (x = 0; x < len1; x++)
			{
				s2 = (*pEdges)[s1][r1][x].getS2();
				r2 = (*pEdges)[s1][r1][x].getR2();
				w2 = (*pEdges)[s1][r1][x].getWeight();
				maxRes[s1][r1] += w2;
				maxRes[s2][r2] += w2;
				pCL->setMaxValue(max(w2, pCL->getMaxValue()));
			}
		}
	}

	for (s1 = 0; s1 < _numSeqs; s1++)
	{
		for (r1 = 1; r1 <= alignPtr->getSeqLength (s1 + 1) && r1 < (*pEdges)[s1].size(); r1++)
		{
			pCL->setMaxExtValue(max(maxRes[s1][r1], pCL->getMaxExtValue()));
		}
	}

	/*
	 * Free maxRes memory
	 */
	for (s1 = 0; s1 < _numSeqs; s1++)
	{
		vector<int>().swap(maxRes[s1]);
	}
	vector<vector<int>>().swap(maxRes);
	
	if (pCL->getNormalise())
	{
		pCL->setNoMatch((pCL->getNoMatch()*pCL->getNormalise())/pCL->getMaxExtValue());
	}
}

void cl2PairListExt (clustalw::Alignment *alignPtr, vector<int> *group)
{
	vector<vector<int>>	ls(2);
	vector<int>		ll(2);
	int 			l1,
				l2,
				len1,
				len2,
				a,
				b,
				g,
				s,
				r,
				i,
				t_s,
				t_r,
				t_w,
				t_s2,
				t_r2,
				t_w2,
				si,
				p1,
				p2,
				s1,
				r1,
				gapPos1,
				gapPos2,
				norm2,
				nused;
	float 			nscore, 
				score, 
				tot, 
				filter, 
				avg=0;
	clustalw::Constraints	*pCL;
	vector<ExtEdge>		*extEdges;
	vector<int>		*indExtEdges;
	vector<vector<vector<clustalw::Edge>>>	*pEdges;
	const int		NORM_F=1000.0;
	const vector<int>	*ptrSeq;

/*
	pos=aln2pos_simple ( A,-1, ns, ls);
	inv_pos=(int**)vcalloc ((CL->S)->nseq, sizeof (int*));
	for (a=0; a<ns[1]; a++)inv_pos[ls[1][a]] =seq2inv_pos(A->seq_al[ls[1][a]]);
*/

	pCL = alignPtr->getPCL();
//	pCL->extEdges = new vector<ExtEdge>();
//	pCL->indExtEdges = new vector<int>;
	pEdges = pCL->getEdges();

	for (int i=1; i < group->size(); i++)
	{
		if ((*group)[i] == 0)
			continue;
		
		if ((*group)[i] == 1)
			ls[0].push_back(i-1);
//			ls[0].push_back((int)alignPtr->getUniqueId(i) - 1);
		else
			ls[1].push_back(i-1);
//			ls[1].push_back((int)alignPtr->getUniqueId(i) - 1);
	}

	ll[0] = l1 = alignPtr->getSeqLength (ls[0][0] + 1);
	ll[1] = l2 = alignPtr->getSeqLength (ls[1][0] + 1);

	gapPos1 = clustalw::userParameters->getGapPos1();
	gapPos2 = clustalw::userParameters->getGapPos2();

	vector<int>	sl1 (alignPtr->getNumSeqs(), -1);
	vector<int>	sl2 (alignPtr->getNumSeqs(), -1);
	vector<float>	norm (l1+1);

//	printf ("\n Profile 1 sequences: ");
	for (a = 0; a < ls[0].size(); a++)
	{
		int	seqId;
		seqId = (int)(alignPtr->getUniqueId(ls[0][a]+1)) - 1;
		sl1[seqId]=a;
//		sl1[ls[0][a]]=a;
//		printf ("%d ", ls[0][a]);
//		printf ("%d ", seqId);
	}

	vector<vector<int>>	pos(ls[0].size(), vector<int> (l1 + 1, 0));

	for (a = 0; a < ls[0].size(); a++)
	{
		ptrSeq = alignPtr->getSequence(ls[0][a] + 1);
//		printf ("\n Seq %d: ", ls[0][a]);
//		printf ("\n Seq %d: ", (int)alignPtr->getUniqueId(ls[0][a]+1) - 1);
		for (g = 1, r = 1; g < ptrSeq->size(); g++)
		{
			b = ptrSeq->at(g);
//			printf ("%02d ", b);
			
			if (b != gapPos1 && b != gapPos2)
				pos[a][g] = r++;
		}
	}
//	printf ("\n");
//
//	printf ("\n Profile 2 sequences: ");
	for (a = 0; a < ls[1].size(); a++)
	{
		int	seqId;
		seqId = (int)(alignPtr->getUniqueId(ls[1][a]+1)) - 1;
		sl2[seqId]=a;
//		sl2[ls[1][a]]=a;
//		printf ("%d ", ls[1][a]);
//		printf ("%d ", seqId);
	}

	vector<vector<int>>	inv_pos(ls[1].size(), vector<int> (l2 + 1, 0));

	for (a = 0; a < ls[1].size(); a++)
	{
		ptrSeq = alignPtr->getSequence(ls[1][a] + 1);
//		printf ("\n Seq %d: ", (int)alignPtr->getUniqueId(ls[1][a]+1) - 1);
//		printf ("\n Seq %d: ", ls[1][a]);
		for (g = 1, r = 1; g < ptrSeq->size(); g++)
		{
			b = ptrSeq->at(g);
//			printf ("%02d ", b);

			if (b != gapPos1 && b != gapPos2)
				inv_pos[a][r++] = g;
		}
	}
//	printf ("\n");

	//data structure used for norm 2

/*
	vector<vector<int>>	nr (2, vector<int>(max(l1, l2)+1, 0));

	for (g=0; g<2; g++)
	{
		for (a=0; a<ll[g]; a++)
		{
			nr[g][a+1]=0;

			for (b = 0; b < ls[g].size(); b++)
				 // Descobrir como identificar gap no profile
				 //if (A->seq_al[ls[g][b]][a]!='-')nr[g][a+1]++;
				nr[g][a+1]++;
		}
	}
*/

	vector<vector<float>>	used (l2 + 1, vector<float>(2, 0));
	vector<int>		used_list (l2 + 1);
	nused = 0;
	norm2 = 0;

	extEdges = new vector<ExtEdge>;
	extEdges->reserve (l1 + 1);

	for (p1 = 1; p1 <= l1; p1++)
	{
		int norm1 = 0;

		for (tot=0,nused=0,si=0; si < ls[0].size(); si++)
		{
//			s = ls[0][si];
			s = (int)(alignPtr->getUniqueId(ls[0][si]+1)) - 1;
			//r = pos[s][p1-1]; **** Tenho que entender esse pos
			//r = p1;
			r = pos[si][p1];

			// pos is a matrix (seqs at group1, len group1) to convert from alignment to original seqs
			// pos[si][p1] == 0  means that there is a gap at p1 of seq si
			if (!r)
				continue;

			len1 = (*pEdges)[s][r].size();

			for (a = 0; a < len1; a++)
			{
				t_s = (*pEdges)[s][r][a].getS2();
				t_r = (*pEdges)[s][r][a].getR2();
				t_w = (*pEdges)[s][r][a].getWeight();

				if (sl1[t_s] != -1)
					continue;//do not extend within a profile

				norm1++;
				norm[p1]++; 
				norm2=0;
				len2 = (*pEdges)[t_s][t_r].size();

				for (b = -1; b < len2; b++)
				{
					if (b == -1)
					{
						t_s2 = t_s;
						t_r2 = t_r;
						t_w2 = t_w;
					}
					else
					{
						t_s2 = (*pEdges)[t_s][t_r][b].getS2();
						t_r2 = (*pEdges)[t_s][t_r][b].getR2();
						t_w2 = (*pEdges)[t_s][t_r][b].getWeight();
					}

					if (sl2[t_s2] != -1)
					{
						// tenho que entender esse inv_pos
						p2 = inv_pos[sl2[t_s2]][t_r2];
						//p2 = t_r2;
						score = min(((float)t_w/(float)NORM_F),((float)t_w2/(float)NORM_F));

						if (!used[p2][1] && score > 0)
						{
							used_list[nused++] = p2;
						}

						norm2++;
						tot += score;
						used[p2][0] += score;
						used[p2][1]++;
					}
				}
			}
		}

		filter=0.01;
		for (a = 0; a < nused; a++)
		{
			p2 = used_list[a];
			nscore = used[p2][0]/tot; //Normalized score used for filtering
			score = used[p2][0];
			used[p2][0] = used[p2][1] = 0;

			if (nscore > filter && p1!=0 && p2!=0 && p1!=l1 && p2!=l2)
			{
/*
				if (normalisation_mode==0)
				{ 
*/
					score=((norm[p1]>0)?score/norm[p1]:0);
/*
				}
				else if (normalisation_mode==1)
				{
					score/=(float)((CL->S)->nseq*nr[0][p1]*nr[1][p2]);
				}
				else if (normalisation_mode==2)
				{
					score/=(CL->S)->nseq;
				}
*/

				score*=NORM_F;
				//printf ("%d %d %d %f ", p1, p2, ((l1-(p1))+(p2)), score);
				try
				{
//					clustalw::ExtEdge	e (p1, p2, ((l1-(p1))+(p2)), score);
					extEdges->push_back ({p1, p2, ((l1-(p1))+(p2)), score});
				}
				catch (exception& ex)
				{
					printf ("Não consegui inserir na lista\n");
					exit (1);
				}
			}
		}
	}
//	printf ("\n");

/*
	for (i = 0; i < extEdges->size(); i++)
		printf ("(%d %d %d %f) ", (*extEdges)[i].p1, (*extEdges)[i].p2, (*extEdges)[i].diag, (*extEdges)[i].score);

	printf ("\n");
*/
	/*
	 * MJ
	 * acho que tenho que liberar o extEdges anterior
	 * O mesmo para o indExtEdges
	 */
	pCL->extEdges = extEdges;
	p1 = (*pCL->extEdges)[0].p1;
	indExtEdges = new vector<int>;
	indExtEdges->reserve (l1 + 1);

	while (indExtEdges->size() < p1-1)
		indExtEdges->push_back(-1);

	indExtEdges->push_back(0);
	if (indExtEdges->size() != p1)
	{
		printf ("Error at indExtEdge!!!\n");
		exit (1);
	}

	for (int i = 1; i < extEdges->size(); i++)
	{
		if ((*extEdges)[i].p1 == p1)
			continue;

		p1 = (*extEdges)[i].p1;
		while (indExtEdges->size() < p1-1)
			indExtEdges->push_back(-1);

		indExtEdges->push_back(i);
//		printf ("(%lu, %d)", indExtEdges->size() -1, i);

		if (indExtEdges->size() != p1)
		{
			printf ("Error at indExtEdge!!!Size = %lu, p1 = %d\n", indExtEdges->size(), p1);
			exit (1);
		}
	}

	pCL->indExtEdges = indExtEdges;
/*
	for (int i = 0; i < pCL->indExtEdges->size(); i++)
		printf ("(%d - %d) ", i, (*pCL->indExtEdges)[i]);
	printf ("\n");
*/

/*
  free_float (used, -1);
  vfree (used_list);
  free_int (inv_pos, -1);
  free_int (pos, -1);
  vfree (sl2);vfree (sl1);
  vfree(norm);
  return n_in[0];
*/
}
