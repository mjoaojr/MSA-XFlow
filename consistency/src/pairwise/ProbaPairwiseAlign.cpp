/**
 * Author: Mario Joao Jr.
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "../consistency/Constraints.h"
#include "../consistency/ProbConsConsistency.h"
#include "../consistency/TCoffeeConsistency.h"
#include "ProbaPairwiseAlign.h"
#include <math.h>
#include <time.h>
#include <algorithm>

namespace clustalw
{

ProbaPairwiseAlign::ProbaPairwiseAlign()
{
	_maxAlnLength = 0;
	pemitPairs = NULL;
	TinsProb_ml = 0; 
	TmatchProb_ml = 0;
}

void ProbaPairwiseAlign::pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, int iEnd, int jStart, int jEnd)
{
	int si, sj, i;
	double _score;
	int _matAvgScore;
	int maxRes;
	int 		_numSeqs = alignPtr->getNumSeqs();
	Constraints	*pCL = new Constraints(_numSeqs);
	clock_t         ini, fim;
    
	try
	{
        
		if(distMat->getSize() != alignPtr->getNumSeqs() + 1)
		{
		    cerr << "The distance matrix is not the right size!\n"
			 << "Need to terminate program.\n";
		    exit(1);
		}
		if((iStart < 0) || (iEnd < iStart) || (jStart < 0) || (jEnd < jStart))
		{
		    cerr << "The range for pairwise Alignment is incorrect.\n"
			 << "Need to terminate program.\n";
		    exit(1);
		}
		
		_maxAlnLength = alignPtr->getMaxAlnLength();
	    
		if(_numSeqs == 0)
		{
		    return;
		}
	    
		int num = (2 * _maxAlnLength) + 1;
	    
		PairScaleValues scaleValues;
		maxRes = subMatrix->getPairwiseMatrix(matrix, scaleValues, _matAvgScore);
		if (maxRes == 0)
		{
		    cerr << "Could not get the substitution matrix\n";
		    return;
		}

		alignPtr->setConstraints(pCL);
		align = alignPtr;

		ini = clock();
		for (si = utilityObject->MAX(0, iStart); si < _numSeqs && si < iEnd; si++)
		{
		    lenSeq1 = alignPtr->getSeqLength(si + 1);

		    for (sj = utilityObject->MAX(si+1, jStart+1); sj < _numSeqs && sj < jEnd; sj++)
		    {
			lenSeq2 = alignPtr->getSeqLength(sj + 1);
			// align the sequences
		    
			_ptrToSeq1 = alignPtr->getSequence(si + 1);
			_ptrToSeq2 = alignPtr->getSequence(sj + 1);
		    
			// Use Proba Pair to calculate score
			probaScore = proba_pair_wise(si, sj);
			distMat->SetAt(si + 1, sj + 1, probaScore);
			distMat->SetAt(sj + 1, si + 1, probaScore);
			
			if(userParameters->getDisplayInfo())
			{
			    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %.4f",
						si+1, sj+1, (float)probaScore);     
			}
		    }
		}
		fim = clock();
		printf ("\nTempo PROBA: %g\n",  ((double) (fim - ini)) / CLOCKS_PER_SEC);

/*
		if (!userParameters->getConsistency())
			return;
*/

		ini = clock();
                relaxCL (align);
                fim = clock();
                printf ("\nTempo relax: %g\n",  ((double) (fim - ini)) / CLOCKS_PER_SEC);
                maxConstraints (align);
//		pCL->printEdges();

		/*
		 * Just one ConsistencyRep
		pCL = align->getPCL();
		vector<vector<SparseMatrix *> > sm = pCL->getSparseMatrices();
		vector<vector<SparseMatrix *> > newSparseMatrices = DoRelaxation (align, sm);

		int numSeqs = align->getNumSeqs();
		// now replace the old posterior matrices
		for (int i = 0; i < _numSeqs; i++)
		{
			for (int j = 0; j < _numSeqs; j++)
			{
				delete sm[i][j];
				sm[i][j] = newSparseMatrices[i][j];
			}
		}
		 */

		/*
		 * Entender como sm se relaciona com o alinhamento... (30/10/2023)
		 * A CL de um par de sequências é uma Sparse Matrix? (01/11/2023)
		 * Acho que não serve para nada (05/11/2023)
		pCL->setSparseMatrices (sm);
		 */

	}
	catch(const exception& e)
	{
		cerr << "An exception has occured in the ProbaPairwiseAlign class.\n"
		     << e.what() << "\n";
		exit(1);    
	}
}

inline int ProbaPairwiseAlign::getElem1(int i)
{
    return (*_ptrToSeq1)[i+1];
}

inline int ProbaPairwiseAlign::getElem2(int i)
{
    return (*_ptrToSeq2)[i+1];
}

float ProbaPairwiseAlign::LOOKUP (float x)
{
	if (x <= 1.00f) return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
	if (x <= 2.50f) return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
	if (x <= 4.50f) return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;

	return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

/*********************************************************************/
/*                                                                   */
/*                         Fast Log Additions (adapted from Probcons)*/
/*                                                                   */
/*                                                                   */
/*********************************************************************/
double EXP (double x)
{
	if (x > -2)
	{
		if (x > -0.5)
		{
			if (x > 0)
				return exp(x);
			return (((0.03254409303190190000*x + 0.16280432765779600000)*x + 0.49929760485974900000)*x + 0.99995149601363700000)*x + 0.99999925508501600000;
		}
		if (x > -1)
			return (((0.01973899026052090000*x + 0.13822379685007000000)*x + 0.48056651562365000000)*x + 0.99326940370383500000)*x + 0.99906756856399500000;
		return (((0.00940528203591384000*x + 0.09414963667859410000)*x + 0.40825793595877300000)*x + 0.93933625499130400000)*x + 0.98369508190545300000;
	}

	if (x > -8)
	{
		if (x > -4)
			return (((0.00217245711583303000*x + 0.03484829428350620000)*x + 0.22118199801337800000)*x + 0.67049462206469500000)*x + 0.83556950223398500000;
		return (((0.00012398771025456900*x + 0.00349155785951272000)*x + 0.03727721426017900000)*x + 0.17974997741536900000)*x + 0.33249299994217400000;
	}

	if (x > -16)
		return (((0.00000051741713416603*x + 0.00002721456879608080)*x + 0.00053418601865636800)*x + 0.00464101989351936000)*x + 0.01507447981459420000;
	return 0;
}


void ProbaPairwiseAlign::LOG_PLUS_EQUALS (float *x, float y)
{
	if (x[0] < y)
		x[0] = (x[0] == LOG_ZERO || y - x[0] >= LOG_UNDERFLOW_THRESHOLD) ? y : LOOKUP(y-x[0]) + x[0];
	else
		x[0] = (y == LOG_ZERO || x[0] - y >= LOG_UNDERFLOW_THRESHOLD) ? x[0]  : LOOKUP(x[0]-y) + y;
}

float ProbaPairwiseAlign::LOG_ADD (float x, float y)
{
	if (x < y) 
		return (x == LOG_ZERO || y - x >= LOG_UNDERFLOW_THRESHOLD) ? y : LOOKUP((y-x)) + x;
	return (y == LOG_ZERO || x - y >= LOG_UNDERFLOW_THRESHOLD) ? x : LOOKUP((x-y)) + y;
}

float **ProbaPairwiseAlign::declare_float (int rowCount, int colCount)
{
	float** a = new float*[rowCount];
	for(int i = 0; i < rowCount; ++i)
		a[i] = new float[colCount];

	return a;
}

float ***declare_array3 (int dim1, int dim2, int dim3, size_t size)
{
	int a1, a2;
	float ***p;

	p=(float***)malloc (dim1 * sizeof (float**));
	for (a1 = 0; a1 < dim1; a1++)
	{
		p[a1] = (float**)malloc (dim2 * sizeof (float*));
		for (a2 = 0; a2 < dim2; a2++)
		{
			p[a1][a2] = (float*)malloc (dim3 * size);
		}
	}
	return p;
}

void free_array3 (float ***p, int dim1, int dim2)
{
	int a1, a2;

	for (a1 = 0; a1 < dim1; a1++)
	{
		for (a2 = 0; a2 < dim1; a2++)
		{
			free (p[a1][a2]);
		}
		free (p[a1]);
	}
	free (p);
}


void ProbaPairwiseAlign::get_emitPairs (const char *alp, vector<float> s)
{
	static char *rmat;
	float k=0, t=0;
	int a, b, c, l;
	int **M;

	l=strlen (alp);

	k=log (2)/2;
	for (a=0; a<l; a++)
	for (b=0; b<l; b++)
	{
		int sc;
		float e;
		e=s[a]*s[b];
		sc=matrix[alp[a]-'A'][alp[b]-'A'];
		pemitPairs[a][b]=e*exp ((double)sc*k);
	}

	for (a=0; a<l; a++)
		for (b=0; b<l; b++)
			t+=pemitPairs[a][b];

	for (a=0; a<l; a++)
		for (b=0; b<l; b++)
			pemitPairs[a][b]=pemitPairs[a][b]/t;

	t=0;

	for (a=0; a<l; a++)
		for (b=0; b<l; b++)
			t+=pemitPairs[a][b];
}

int ProbabilisticModel (int NumMatrixTypes, int NumInsertStates, const float *initDistribMat, vector<float> emitSingle, float **emitPairs, const float *gapOpen, const float *gapExtend, float **transMat, vector<float> &initialDistribution, float **matchProb, float **insProb, float **transProb)
{
	// build transition matrix
	int i, j;

	//Maybe an Issue with this topology

	transMat[0][0] = 1;
	for (i = 0; i < NumInsertStates; i++)
	{
		transMat[0][2*i+1] = gapOpen[2*i];
		transMat[0][2*i+2] = gapOpen[2*i+1];
		transMat[0][0] -= (gapOpen[2*i] + gapOpen[2*i+1]);

		transMat[2*i+1][2*i+1] = gapExtend[2*i];
		transMat[2*i+2][2*i+2] = gapExtend[2*i+1];
		transMat[2*i+1][2*i+2] = 0;
		transMat[2*i+2][2*i+1] = 0;
		transMat[2*i+1][0] = 1 - gapExtend[2*i];
		transMat[2*i+2][0] = 1 - gapExtend[2*i+1];
	}

	// create initial and transition probability matrices
	for (i = 0; i < NumMatrixTypes; i++)
	{
		initialDistribution[i] = (float)log ((float)initDistribMat[i]);
		for (j = 0; j < NumMatrixTypes; j++)
			transProb[i][j] = (float)log ((float)transMat[i][j]);
	}

	// create insertion and match probability matrices
	for (i = 0; i < 256; i++)
	{
		for (j = 0; j < NumMatrixTypes; j++)
		{
			insProb[i][j] = (float)log((float)emitSingle[i]);
		}
		for (j = 0; j < 256; j++)
		{
			matchProb[i][j] = (float)log((float)emitPairs[i][j]);
		}
	}

/* DEBUG MJ
	for (i = 'A'; i <= 'Z'; i++)
	{
		for (j = 'A'; j < i; j++)
		{
			printf ("%.7f ", matchProb[i][j]);
		}
		printf ("\n");
	}
	exit (0);
*/

	return 1;
}

int ProbaPairwiseAlign::get_tot_prob2 (int nstates, float **matchProb, float **insProb, float *TmatchProb, float ***TinsProb, int mode)
{
	static double **prf1;
	static double **prf2;
	int i, j, k, ij,r,r1,r2; 
 
/*
 * Não vi serventia
	if (mode==SEQUENCE)
	{
		int s1, s2, a;
		int *nns, **nls;
		Alignment *NA1, *NA2;
		char *sst1;
		char *sst2;

		nns=(int*)vcalloc ( 2, sizeof (int));
		nls=(int**)vcalloc (2, sizeof (int*));

		s1=A1->order[ls[0][0]][0];
		s2=A2->order[ls[1][0]][0];
		NA1=seq2R_template_profile (CL->S,s1);
		NA2=seq2R_template_profile (CL->S,s2);

		sst1=seq2T_template_string((CL->S),s1);
		sst2=seq2T_template_string((CL->S),s2);

		if (NA1 || NA2)
		{
			if (NA1)
			{
				nns[0]=NA1->nseq;
				nls[0]=(int*)vcalloc (NA1->nseq, sizeof (int));
				for (a=0; a<NA1->nseq; a++)
					nls[0][a]=a;
				NA1->seq_al[NA1->nseq]=sst1;
				sprintf (NA1->name[NA1->nseq], "sst1");
			}
			else
			{
				NA1=A1;
				nns[0]=ns[0];
				nls[0]=(int*)vcalloc (ns[0], sizeof (int));
				for (a=0; a<ns[0]; a++)
					nls[0][a]=ls[0][a];
			}

			if (NA2)
			{
				nns[1]=NA2->nseq;
				nls[1]=(int*)vcalloc (NA2->nseq, sizeof (int));
				for (a=0; a<NA2->nseq; a++)
					nls[1][a]=a;
				NA2->seq_al[NA2->nseq]=sst2;
				sprintf (NA2->name[NA2->nseq], "sst2");
			}
			else
			{
				NA2=A2;
				nns[1]=ns[1];
				nls[1]=(int*)vcalloc (ns[1], sizeof (int));
				for (a=0; a<ns[1]; a++)
					nls[1][a]=ls[1][a];
			}

			get_tot_prob2 (NA1, NA2, nns, nls, nstates, matchProb, insProb, TmatchProb, TinsProb, CL,PROFILE);
			vfree (nns); 
			free_int (nls,-1);
			return 1;
		}
	}
	else if (mode == PROFILE)
	{
		static int display_mode;
		if (!display_mode)
		{
			fprintf ( stderr, "\n! profile/profile alignment --- use proba_pair\n");
			display_mode=1;
		}
	}
*/
        
	/*
	 *	Constante na propria classe
	lu=dirichlet_code2aa_lu();
	 */
  
	/*
	 * Não preciso disso
	prf1=aln2prf (A1, ns[0], ls[0], lenSeq1, prf1);
	prf2=aln2prf (A2, ns[1], ls[1], lenSeq2, prf2);
	*/

    
	//get Ins for I
	for (i=1; i<=lenSeq1; i++)
	{
//              printf ("%c(%d) - ", lu[convlu[getElem1(i-1)]], getElem1(i-1));
		for (k=0; k<nstates; k++)
		{
			TinsProb[0][k][i] = (float)insProb[asclu[getElem1(i-1)]][k];
		}
//              printf ("\n");
	}

	//get Ins for J
	for (j=1; j<=lenSeq2; j++)
	{
		for (k=0; k<nstates; k++)
		{
			TinsProb[1][k][j] = (float)insProb[asclu[getElem2(j-1)]][k];
		}
	}

	for (ij=0,i=0; i<=lenSeq1; i++)
		for (j=0; j<=lenSeq2; j++, ij++)
		{
			float tot=0,f;

			if (i==0 || j==0)
				continue;

			TmatchProb[ij]=0;
			for (tot=0,r1=0; r1<22; r1++)
			{
				for (r2=0; r2<22; r2++)
				{
					/*
					f=(float)prf1[i-1][r1]*(float)prf2[j-1][r2];
					TmatchProb[ij]+=matchProb[lu[r1]][lu[r2]]*f;
					*/
					if ((lu[convlu[getElem1(i-1)]] == lu[r1]) && (lu[convlu[getElem2(j-1)]] == lu[r2]))
					{
						TmatchProb[ij]+=matchProb[lu[r1]][lu[r2]];
//						printf ("i=%02d,j=%02d,r1=%02d,r2=%02d(%d)=%.7f[%.7f]\n", i, j, r1, r2, ij, TmatchProb[ij],matchProb[lu[r1]][lu[r2]]);
					}
				}
			}
		}

/*
	for (ij=0,i=0; i<=lenSeq1; i++)
	{
		for (j=0; j<=lenSeq2; j++, ij++)
		{
			if (i==0 || j==0)
				continue;

			TmatchProb[ij]+=matchProb[asclu[getElem1(i-1)]][asclu[getElem2(j-1)]];
		}
	}

	for (i = 0; i < ij; i++)
		printf ("%.7f ", TmatchProb[ij]);
	printf ("\n");
*/
  
  return 1;
}

float * ProbaPairwiseAlign::forward_proba_pair_wise (int NumMatrixTypes, int NumInsertStates, float **transMat, vector<float> initialDistribution,float *matchProb, float ***insProb, float **transProb)
{
	float *forward;
	int max_l;
	int k, i, j,ij, i1j1, i1j, ij1, m;
	int l,a;

	l=(lenSeq1+1)*(lenSeq2+1)*NumMatrixTypes;

	forward = (float *)malloc (l * sizeof(float));
	max_l=l;

	for (a=0; a<l; a++)forward[a]=LOG_ZERO;

	forward[0 + NumMatrixTypes * (1 * (lenSeq2+1) + 1)] = initialDistribution[0] + matchProb[lenSeq2+2];

	for (k = 0; k < NumInsertStates; k++)
	{
		forward[2*k+1 + NumMatrixTypes * (1 * (lenSeq2+1) + 0)] = initialDistribution[2*k+1] + insProb[0][k][1];
		forward[2*k+2 + NumMatrixTypes * (0 * (lenSeq2+1) + 1)] = initialDistribution[2*k+2] + insProb[1][k][1];
	}

	// remember offset for each index combination
	ij = 0;
	i1j = -lenSeq2 - 1;
	ij1 = -1;
	i1j1 = -lenSeq2 - 2;

	ij *= NumMatrixTypes;
	i1j *= NumMatrixTypes;
	ij1 *= NumMatrixTypes;
	i1j1 *= NumMatrixTypes;


	// compute forward scores
	for (m=0,i = 0; i <= lenSeq1; i++)
	{
		for (j = 0; j <= lenSeq2; j++, m++)
		{
			if (i > 1 || j > 1)
			{
				if (i > 0 && j > 0)
				{
					//Sum over all possible alignments
					forward[0 + ij] = forward[0 + i1j1] + transProb[0][0];
					for (k = 1; k < NumMatrixTypes; k++)
					{
						LOG_PLUS_EQUALS (&forward[0 + ij], forward[k + i1j1] + transProb[k][0]);
					}
					forward[0 + ij] += matchProb[m];
				}
				if ( i > 0)
				{
					for (k = 0; k < NumInsertStates; k++)
					{
						forward[2*k+1 + ij] = insProb[0][k][i] + LOG_ADD (forward[0 + i1j] + transProb[0][2*k+1],forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1]);
					}
				}
				if (j > 0)
				{
					for (k = 0; k < NumInsertStates; k++)
					{
						forward[2*k+2 + ij] = insProb[1][k][j] +LOG_ADD (forward[0 + ij1] + transProb[0][2*k+2],forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2]);
					}
				}
			}

			ij += NumMatrixTypes;
			i1j += NumMatrixTypes;
			ij1 += NumMatrixTypes;
			i1j1 += NumMatrixTypes;
		}
	}

/*
	for (i = 0; i < l; i+=NumMatrixTypes)
	{
		printf ("%f %f %f %f %f\n", forward[i], forward[i+1], forward[i+2], forward[i+3], forward[i+4]);
	}
*/
	return forward;
}

float *ProbaPairwiseAlign::backward_proba_pair_wise (int NumMatrixTypes, int NumInsertStates, float **transMat, vector<float> initialDistribution,float *matchProb, float ***insProb, float **transProb)
{
	float *backward;
	int max_l;

	int k, i, j,ij, i1j1, i1j, ij1,a, l, m;

	l=(lenSeq1+1)*(lenSeq2+1)*NumMatrixTypes;

	backward = (float *)malloc (l * sizeof(float));
	max_l=l;

	for (a=0; a<l; a++)backward[a]=LOG_ZERO;

	for (k = 0; k < NumMatrixTypes; k++)
		backward[NumMatrixTypes * ((lenSeq1+1) * (lenSeq2+1) - 1) + k] = initialDistribution[k];
  
	//Difference with Probcons: this emission is not added to the bward
	backward[NumMatrixTypes * ((lenSeq1+1) * (lenSeq2+1) - 1) + 0]+=matchProb[(lenSeq1+1) * (lenSeq2+1) - 1];
	// remember offset for each index combination
	ij = (lenSeq1+1) * (lenSeq2+1) - 1;

	i1j = ij + lenSeq2 + 1;
	ij1 = ij + 1;
	i1j1 = ij + lenSeq2 + 2;
	ij *= NumMatrixTypes;
	i1j *= NumMatrixTypes;
	ij1 *= NumMatrixTypes;
	i1j1 *= NumMatrixTypes;

	// compute backward scores
	for (i = lenSeq1; i >= 0; i--)
	{
		for (j = lenSeq2; j >= 0; j--)
		{
			if (i < lenSeq1 && j < lenSeq2)
			{
				m=((i+1)*(lenSeq2+1))+j+1;//The backward and the forward are offset by 1
				float ProbXY = backward[0 + i1j1] + matchProb[m];

				for (k = 0; k < NumMatrixTypes; k++)
				{
					LOG_PLUS_EQUALS (&backward[k + ij], ProbXY + transProb[k][0]);
				}
			}
			if (i < lenSeq1)
			{
				for (k = 0; k < NumInsertStates; k++)
				{
					LOG_PLUS_EQUALS (&backward[0 + ij], backward[2*k+1 + i1j] + insProb[0][k][i+1] + transProb[0][2*k+1]);
					LOG_PLUS_EQUALS (&backward[2*k+1 + ij], backward[2*k+1 + i1j] + insProb[0][k][i+1] + transProb[2*k+1][2*k+1]);
				}
			}
			if (j < lenSeq2)
			{
				for (k = 0; k < NumInsertStates; k++)
				{
					//+1 because the backward and the forward are offset by 1
					LOG_PLUS_EQUALS (&backward[0 + ij], backward[2*k+2 + ij1] + insProb[1][k][j+1] + transProb[0][2*k+2]);
					LOG_PLUS_EQUALS (&backward[2*k+2 + ij], backward[2*k+2 + ij1] + insProb[1][k][j+1] + transProb[2*k+2][2*k+2]);
				}
			}

			ij -= NumMatrixTypes;
			i1j -= NumMatrixTypes;
			ij1 -= NumMatrixTypes;
			i1j1 -= NumMatrixTypes;
		}
	}
 
	return backward;
}

float ProbaPairwiseAlign::ComputeTotalProbability (int seq1Length, int seq2Length,int NumMatrixTypes, int NumInsertStates,float *forward, float *backward)
{
	float totalForwardProb = LOG_ZERO;
	float totalBackwardProb = LOG_ZERO;
	int k;

	for (k = 0; k < NumMatrixTypes; k++)
	{
		LOG_PLUS_EQUALS (&totalForwardProb,forward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] + backward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
//		printf ("f[%d]=%f b[%d]=%f\n", k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1), forward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)], k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1), backward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
	}

	totalBackwardProb =forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] +backward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)];

	for (k = 0; k < NumInsertStates; k++)
	{
		LOG_PLUS_EQUALS (&totalBackwardProb,forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] +backward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)]);
		LOG_PLUS_EQUALS (&totalBackwardProb,forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] +backward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)]);
	}

//	printf ("%f %f\n", totalForwardProb, totalBackwardProb);
    
	return (totalForwardProb + totalBackwardProb) / 2;
}

bool customCmp(vector<int> a, vector<int> b) { return a[2] > b[2]; }

//Constraint_list *ProbaMatrix2CL (Alignment *A, int *ns, int **ls, int NumMatrixTypes, int NumInsertStates, float *forward, float *backward, float thr, Constraint_list *CL)
float ProbaPairwiseAlign::ProbaMatrix2CL (int NumMatrixTypes, int NumInsertStates, float *forward, float *backward, float thr, int s1, int s2, vector<float> &posterior)
{
	float totalProb,
		maxTotProb;
	int ij, i, j,k;
	float *twoRows = new float[(lenSeq2+2)*2];
	float *oldRow = twoRows;
	float *newRow = twoRows + lenSeq2 + 2;
	vector <vector<int>>	list;
	static int list_max;
	int sim;
	int list_size;
	int list_n;
	int old_n=0;
	double v;
	int a;
	static float F=4; //potential number of full suboptimal alignmnents incorporated in the library
	static int tot_old, tot_new;
	Constraints	*pCL;
  
	list_size=lenSeq1*lenSeq2;

	totalProb = ComputeTotalProbability (lenSeq1,lenSeq2,NumMatrixTypes, NumInsertStates,forward, backward);

	vector<float> *posteriorPtr = new vector<float>((lenSeq1+1) * (lenSeq2+1));
	posterior = *posteriorPtr;
	vector<float>::iterator ptr = posterior.begin();

	// initialization
	for (i = 0; i <= lenSeq2+1; i++)
	{
		oldRow[i] = 0;
	}

	pCL = align->getConstraints();
	ij = 0;
	for (list_n=0,ij=0,i =0; i <= lenSeq1; i++)
	{
		// initialize left column
		newRow[0] = 0;
 
		for (j =0; j <= lenSeq2; j++, ij+=NumMatrixTypes)
		{
			v= EXP (min(LOG_ONE,(forward[ij] + backward[ij] - totalProb)));
			*(ptr++) = v;
			newRow[j+1] = ChooseBestOfThree (v + oldRow[j], newRow[j], oldRow[j+1]);
			if (v>thr)//Conservative reduction of the list size to speed up the sorting
			{
/* --
				vector<int> vet;
				vet.push_back (i);
				vet.push_back (j);
				vet.push_back ((int)((float)v*(float)NORM_F));
				list.push_back (vet);
				list_n++;
//				printf ("(%d, %d, %d) ", vet[0], vet[1], vet[2]);
*/

				ConstrEntry	entry(s1, s2, i, j, (int)((float)v*(float)NORM_F), 1, -1);
				pCL->addEntry(entry);
				list_n++;
			}

			if (v>0.01)
				old_n++;
		}

//		printf ("%d - %f\n", i, newRow[lenSeq2+1]);

		// swap rows
		float *temp = oldRow;
		oldRow = newRow;
		newRow = temp;
	}
//	printf ("\nlist_n = %d\n", list_n);

	posterior[0] = 0;

	maxTotProb = oldRow[lenSeq2+1] / min (lenSeq1, lenSeq2);
//	printf ("chute = %f, totalProb = %f, list_n=%d, old_n=%d\n", maxTotProb, totalProb, list_n, old_n);

// --	sort (list.begin(), list.end(), customCmp);
//	sort_int_inv (list, 3, 2, 0, list_n-1);
//	if (!entry)entry=(int*)vcalloc ( CL->entry_len+1, CL->el_size);

	list_n=min(list_n,(int)(F*min(lenSeq1,lenSeq2)));

/* --
	for (i=0; i<list_n; i++)
	{
		ConstrEntry	entry(s1, s2, list[i][0], list[i][1], list[i][2], 1, -1);
		entry[SEQ1]=s1;
		entry[SEQ2]=s2;
		entry[R1]  =list[i][0];
		entry[R2]  =list[i][1];
		entry[WE]  =list[i][2];
		entry[CONS]=1;
//		printf ("(%d, %d, %d) ", list[i][0], list[i][1], list[i][2]);
		pCL->addEntry(entry);
	}
//	printf ("\nlist_n = %d\n", list_n);
*/

/*
	pCL->printEdges();
	pCL->print();
	pCL->printResIndex();
 */

	tot_new+=list_n;
	tot_old+=old_n;
	// HERE ("LIB_SIZE NEW: %d (new) %d (old) [%.2f]", list_n, old_n, (float)tot_new/(float)tot_old);
	delete (twoRows);
	delete (posteriorPtr);
	return maxTotProb;
}


//int proba_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
float ProbaPairwiseAlign::proba_pair_wise (int s1, int s2)
{
	static int NumMatrixTypes;
	static int NumInsertStates;
	static float **transMat, **insProb, **matchProb, **transProb, **emitPairs, ***TinsProb;
	float *TmatchProb;
	static vector<float>	emitSingle(256, 0),
				initialDistribution(11, 0);
	float	score;
	vector<float>	posterior;
	Constraints	*pCL;

	float thr=0.01;//ProbCons Default
	int i, j;
	float *F, *B;
	int l;
	const char *alphabet;

	if (!transMat)
	{
		static vector<float> s(strlen(alphabetDefault), 0);
		NumInsertStates=2;
		NumMatrixTypes=5;
/*
		if (atoigetenv ("NOBIPHASIC"))
		{
			NumInsertStates=1;
			NumMatrixTypes=3;
		}
*/
		if (pemitPairs == NULL)
		{
			int l,a,b;
			l=strlen (alphabetDefault);
			pemitPairs=declare_float (l,l);
//			s=(float*)vcalloc (l, sizeof (float));
			for (a=0; a<l; a++)
			{
				s[a]=emitSingleDefault[a];
				for (b=0; b<l; b++)
					pemitPairs[a][b]=emitPairsDefault[a][b];
			}
			
/* DEBUG MJ
			for (a=0; a<l; a++)
			{
				for (b=0; b<l; b++)
					printf ("(%d,%d)%.7f ", a, b, pemitPairs[a][b]);
				printf ("\n");
			}
			exit(0);
*/
		}
// This way, it uses the default emitPairs matrix
//		get_emitPairs (alphabetDefault,s);
		alphabet=alphabetDefault;
		emitPairs=declare_float (256, 256);
//		emitSingle=(float*)vcalloc (256, sizeof (float));
		for (i=0; i<256; i++)
		{
			//emitSingle[i]=1e-5;
			emitSingle[i]=1;
			for (j=0; j<256; j++)
				//emitPairs[i][j]=1e-10;
				emitPairs[i][j]=1;
		}
		l=strlen (alphabet);

		for (i=0; i<l; i++)
		{
			int C1,c1, C2,c2;
			c1=tolower(alphabet[i]);
			C1=toupper(alphabet[i]);
			emitSingle[c1]=s[i];
			emitSingle[C1]=s[i];
			for (j=0; j<=i; j++)
			{
				c2=tolower(alphabet[j]);
				C2=toupper(alphabet[j]);

				emitPairs[c1][c2]=pemitPairs[i][j];
				emitPairs[C1][c2]=pemitPairs[i][j];
				emitPairs[C1][C2]=pemitPairs[i][j];
				emitPairs[c1][C2]=pemitPairs[i][j];
				emitPairs[c2][c1]=pemitPairs[i][j];
				emitPairs[C2][c1]=pemitPairs[i][j];
				emitPairs[C2][C1]=pemitPairs[i][j];
				emitPairs[c2][C1]=pemitPairs[i][j];
//			printf ("(%d,%d)%.7f ", c1, c2, emitPairs[c1][c2]);
			}
		}
/* DEBUG
			printf ("\n");
			exit (0);
*/


		transMat=declare_float (2*NumInsertStates+1, 2*NumInsertStates+1);
		transProb=declare_float (2*NumInsertStates+1,2* NumInsertStates+1);
		insProb=declare_float (256,NumMatrixTypes);
		matchProb=declare_float (256, 256);
		
		initialDistribution.resize(2*NumMatrixTypes+1);
		ProbabilisticModel (NumMatrixTypes,NumInsertStates,initDistrib2Default, emitSingle,emitPairs,gapOpen2Default,gapExtend2Default, transMat,initialDistribution,matchProb, insProb,transProb);
	}

	l=(lenSeq1+1)*(lenSeq2+1);
	if (l>TmatchProb_ml)
	{
		TmatchProb_ml=l;
	}

	TmatchProb = (float *)malloc (TmatchProb_ml * sizeof (float));

	l=max(lenSeq1,lenSeq2)+1;
	if ( l>TinsProb_ml)
	{
		TinsProb_ml=l;
		if (TinsProb)
			free_array3 (TinsProb, 2, NumMatrixTypes);
		TinsProb=declare_array3 (2, NumMatrixTypes, TinsProb_ml, sizeof (float));
	}

	get_tot_prob2 (NumMatrixTypes, matchProb, insProb, TmatchProb, TinsProb, SEQUENCE);

	F=forward_proba_pair_wise (NumMatrixTypes,NumInsertStates,transMat, initialDistribution,TmatchProb,TinsProb, transProb);
	B=backward_proba_pair_wise (NumMatrixTypes,NumInsertStates,transMat, initialDistribution,TmatchProb,TinsProb, transProb);
	score = ProbaMatrix2CL(NumMatrixTypes,NumInsertStates, F, B, thr, s1, s2, posterior);

	pCL = align->getPCL();
	pCL->setSparseMatrix (s1, s2, lenSeq1, lenSeq2, posterior);


//	printf ("%f\n", score);

	free (TmatchProb);
	free (F);
	free (B);
	//free_proba_pair_wise();
	return score;
}

}
