/**
 * Author: Mario Joao Jr.
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <cassert>

#include "KmersPairwiseAlign.h"

namespace clustalw
{

KmersPairwiseAlign::KmersPairwiseAlign()
{
	const vector<char> nucleotides = {'A', 'C', 'G', 'T'};
        bool _DNAFlag = userParameters->getDNAFlag();

	_maxAlnLength = 0;

	if (_DNAFlag)
	{
		const vector<enum NX> letters = {NX_A, NX_C, NX_G, NX_T};

		for (int i = 0; i < letters.size(); i++)
		{
			CharToLetterEx[nucleotides[i]] = (unsigned char) letters[i];
			CharToLetterEx[(unsigned char)tolower(nucleotides[i])] = (unsigned char) letters[i];
		}

		ResidueGroup = DNAResidueGroup;
		uResidueGroupCount = uDNAResidueGroupCount;
	}
	else
	{
		const vector<enum AX> letters = {AX_A, AX_B, AX_C, AX_D, AX_E, AX_F, AX_G, AX_H, AX_I, AX_K, AX_L, AX_M, AX_N, AX_A, AX_P, AX_Q, AX_R, AX_S, AX_T, AX_A, AX_V, AX_W, AX_X, AX_Y, AX_Z};

		for (int i = 0; i < letters.size(); i++)
		{
			CharToLetterEx[i] = (unsigned char) letters[i];
		}

		ResidueGroup = ProtResidueGroup;
		uResidueGroupCount = uProtResidueGroupCount;
	}
}

void KmersPairwiseAlign::pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, int iEnd, int jStart, int jEnd)
{
    int si, sj, i;
    int n, m, len1, len2;
    double _score;
    
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
    
        int _numSeqs = alignPtr->getNumSeqs();
        if(_numSeqs == 0)
        {
            return;
        }
    
        int num = (2 * _maxAlnLength) + 1;
        bool _DNAFlag = userParameters->getDNAFlag();

	// Convert to letters
	unsigned **Letters = new unsigned *[_numSeqs];
	for (unsigned si = 0; si < _numSeqs; ++si)
	{
                _ptrToSeq1 = alignPtr->getSequence(si + 1);
		const unsigned uSeqLength = alignPtr->getSeqLength(si + 1);
		unsigned *L = new unsigned[uSeqLength];
		Letters[si] = L;
		for (unsigned n = 0; n < uSeqLength; ++n)
		{
			unsigned char c = getElem1(n);
			L[n] = CharToLetterEx[c];
			if (_DNAFlag && (L[n] >= 4))
				L[n] = 4;
		}
	}

	unsigned **uCommonTupleCount = new unsigned *[_numSeqs];
	for (unsigned n = 0; n < _numSeqs; ++n)
	{
		uCommonTupleCount[n] = new unsigned[_numSeqs];
		memset(uCommonTupleCount[n], 0, _numSeqs*sizeof(unsigned));
	}

	const unsigned uPairCount = (_numSeqs*(_numSeqs + 1))/2;
	unsigned uCount = 0;

        for (si = utilityObject->MAX(0, iStart); si < _numSeqs && si < iEnd; si++)
	{
		const unsigned uSeqLength1 = alignPtr->getSeqLength(si + 1);
		if (uSeqLength1 < 5)
			continue;

		const unsigned uTupleCount = uSeqLength1 - 5;
		const unsigned *L = Letters[si];
		CountTuples(L, uTupleCount, Count1);

                if(userParameters->getDisplayInfo())
                {
                    utilityObject->info("K-mer dist pass 1");     
                }

		for (sj = utilityObject->MAX(0, jStart); sj <= si; sj++)
		{
			const unsigned uSeqLength2 = alignPtr->getSeqLength(sj + 1);
			if (uSeqLength2 < 5)
			{
				if (si == sj)
				{
					distMat->SetAt(si + 1, sj + 1, 0.0);
					distMat->SetAt(sj + 1, si + 1, 0.0);
				}
				else
				{
					distMat->SetAt(si + 1, sj + 1, 1.0);
					distMat->SetAt(sj + 1, si + 1, 1.0);
				}
				continue;
			}

			// First pass through seq 2 to count tuples
			const unsigned uTupleCount = uSeqLength2 - 5;
			const unsigned *L = Letters[sj];
			CountTuples(L, uTupleCount, Count2);

			// Second pass to accumulate sum of shared tuples
			// MAFFT defines this as the sum over unique tuples
			// in seq2 of the minimum of the number of tuples found
			// in the two sequences.
			unsigned uSum = 0;
			for (unsigned n = 0; n < uTupleCount; ++n)
			{
				const unsigned uTuple = GetTuple(L, n);
				uSum += min(Count1[uTuple], Count2[uTuple]);

				// This is a hack to make sure each unique tuple counted only once.
				Count2[uTuple] = 0;
			}
			uCommonTupleCount[si][sj] = uSum;
			uCommonTupleCount[sj][si] = uSum;
		}
	}

	uCount = 0;
	if(userParameters->getDisplayInfo())
	{
	    utilityObject->info("K-mer dist pass 2");     
	}

        for (si = utilityObject->MAX(0, iStart); si < _numSeqs && si < iEnd; si++)
	{
                _ptrToSeq1 = alignPtr->getSequence(si + 1);
/*
		Seq &s1 = *(v[uSeq1]);
		const char *pName1 = s1.GetName();
*/

		double dCommonTupleCount11 = uCommonTupleCount[si][si];
		if (0 == dCommonTupleCount11)
			dCommonTupleCount11 = 1;

//		printf ("dCommonTupleCount11 = %.2f\n", dCommonTupleCount11);

		distMat->SetAt(si + 1, si + 1, 0.0);
//		DF.SetDist(uSeq1, uSeq1, 0);
		for (sj = utilityObject->MAX(si+1, jStart+1); sj < _numSeqs && sj < jEnd; sj++)
//		for (unsigned sj = 0; sj < si; ++sj)
		{
			double dCommonTupleCount22 = uCommonTupleCount[sj][sj];
			if (0 == dCommonTupleCount22)
				dCommonTupleCount22 = 1;

//			printf ("dCommonTupleCount22 = %.2f\n", dCommonTupleCount22);
//			printf ("uCommonTupleCount[si][sj] = %u\n", uCommonTupleCount[si][sj]);

//			const double dDist1 = 3.0*(dCommonTupleCount11 - uCommonTupleCount[si][sj]) /dCommonTupleCount11;
//			const double dDist2 = 3.0*(dCommonTupleCount22 - uCommonTupleCount[si][sj]) /dCommonTupleCount22;
			const double dDist1 = (dCommonTupleCount11 - uCommonTupleCount[si][sj]) /dCommonTupleCount11;
			const double dDist2 = (dCommonTupleCount22 - uCommonTupleCount[si][sj]) /dCommonTupleCount22;

//			printf ("dDist1 = %.2f\n", dDist1);
//			printf ("dDist2 = %.2f\n", dDist2);

			// dMinDist is the value used for tree-building in MAFFT
			const double dMinDist = min(dDist1, dDist2);
//			printf ("dMinDist = %.2f\n", dMinDist);
			distMat->SetAt(si + 1, sj + 1, dMinDist);
			distMat->SetAt(sj + 1, si + 1, dMinDist);
//			DF.SetDist(uSeq1, uSeq2, (float) dMinDist);

			if(userParameters->getDisplayInfo())
			{
			    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %.4f", 
						si+1, sj+1, (double)dMinDist);     
			}
			//const double dEstimatedPctId = TupleDistToEstimatedPctId(dMinDist);
			//g_dfPwId.SetDist(uSeq1, uSeq2, dEstimatedPctId);
			// **** TODO **** why does this make score slightly worse??
			//const double dKimuraDist = KimuraDist(dEstimatedPctId);
			//DF.SetDist(uSeq1, uSeq2, dKimuraDist);
		}
	}

	for (unsigned n = 0; n < _numSeqs; ++n)
	{
		delete[] uCommonTupleCount[n];
		delete[] Letters[n];
	}
	delete[] uCommonTupleCount;
	delete[] Letters;
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the LcsPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}

unsigned KmersPairwiseAlign::GetTuple(const unsigned uLetters[], unsigned n)
{
	/*
	assert(uLetters[n] < uResidueGroupCount);
	assert(uLetters[n+1] < uResidueGroupCount);
	assert(uLetters[n+2] < uResidueGroupCount);
	assert(uLetters[n+3] < uResidueGroupCount);
	assert(uLetters[n+4] < uResidueGroupCount);
//	assert(uLetters[n+5] < uResidueGroupCount);
	*/

	unsigned u1 = ResidueGroup[uLetters[n]];
	unsigned u2 = ResidueGroup[uLetters[n+1]];
	unsigned u3 = ResidueGroup[uLetters[n+2]];
	unsigned u4 = ResidueGroup[uLetters[n+3]];
	unsigned u5 = ResidueGroup[uLetters[n+4]];
	unsigned u6 = ResidueGroup[uLetters[n+5]];

	return u6 + u5*6 + u4*6*6 + u3*6*6*6 + u2*6*6*6*6 + u1*6*6*6*6*6;
}

void KmersPairwiseAlign::CountTuples(const unsigned L[], unsigned uTupleCount, unsigned char Count[])
{
	memset(Count, 0, TUPLE_COUNT*sizeof(unsigned char));
	for (unsigned n = 0; n < uTupleCount; ++n)
	{
		const unsigned uTuple = GetTuple(L, n);
		++(Count[uTuple]);
	}
}

/*
void DistKmer4_6(const SeqVect &v, DistFunc &DF)
{
	const unsigned uSeqCount = v.Length();

	// MJ É preciso??????????????
	// Initialize distance matrix to zero
	/*
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			DF.SetDist(uSeq1, uSeq2, 0);
		}

	const unsigned uPairCount = (uSeqCount*(uSeqCount + 1))/2;
	unsigned uCount = 0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &seq1 = *(v[uSeq1]);
		const unsigned uSeqLength1 = seq1.Length();
		if (uSeqLength1 < 5)
			continue;

		const unsigned uTupleCount = uSeqLength1 - 5;
		const unsigned *L = Letters[uSeq1];
		CountTuples(L, uTupleCount, Count1);
#if	TRACE
		{
		Log("Seq1=%d\n", uSeq1);
		Log("Groups:\n");
		for (unsigned n = 0; n < uSeqLength1; ++n)
			Log("%u", ResidueGroup[L[n]]);
		Log("\n");

		Log("Tuples:\n");
		ListCount(Count1);
		}
#endif

		SetProgressDesc("K-mer dist pass 1");
		for (unsigned uSeq2 = 0; uSeq2 <= uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;
			Seq &seq2 = *(v[uSeq2]);
			const unsigned uSeqLength2 = seq2.Length();
			if (uSeqLength2 < 5)
				{
				if (uSeq1 == uSeq2)
					DF.SetDist(uSeq1, uSeq2, 0);
				else
					DF.SetDist(uSeq1, uSeq2, 1);
				continue;
				}

		// First pass through seq 2 to count tuples
			const unsigned uTupleCount = uSeqLength2 - 5;
			const unsigned *L = Letters[uSeq2];
			CountTuples(L, uTupleCount, Count2);
#if	TRACE
			Log("Seq2=%d Counts=\n", uSeq2);
			ListCount(Count2);
#endif

		// Second pass to accumulate sum of shared tuples
		// MAFFT defines this as the sum over unique tuples
		// in seq2 of the minimum of the number of tuples found
		// in the two sequences.
			unsigned uSum = 0;
			for (unsigned n = 0; n < uTupleCount; ++n)
				{
				const unsigned uTuple = GetTuple(L, n);
				uSum += MIN(Count1[uTuple], Count2[uTuple]);

			// This is a hack to make sure each unique tuple counted only once.
				Count2[uTuple] = 0;
				}
#if	TRACE
			{
			Seq &s1 = *(v[uSeq1]);
			Seq &s2 = *(v[uSeq2]);
			const char *pName1 = s1.GetName();
			const char *pName2 = s2.GetName();
			Log("Common count %s(%d) - %s(%d) =%u\n",
			  pName1, uSeq1, pName2, uSeq2, uSum);
			}
#endif
			uCommonTupleCount[uSeq1][uSeq2] = uSum;
			uCommonTupleCount[uSeq2][uSeq1] = uSum;
			}
		}
	ProgressStepsDone();

	uCount = 0;
	SetProgressDesc("K-mer dist pass 2");
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &s1 = *(v[uSeq1]);
		const char *pName1 = s1.GetName();

		double dCommonTupleCount11 = uCommonTupleCount[uSeq1][uSeq1];
		if (0 == dCommonTupleCount11)
			dCommonTupleCount11 = 1;

		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;

			double dCommonTupleCount22 = uCommonTupleCount[uSeq2][uSeq2];
			if (0 == dCommonTupleCount22)
				dCommonTupleCount22 = 1;

			const double dDist1 = 3.0*(dCommonTupleCount11 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount11;
			const double dDist2 = 3.0*(dCommonTupleCount22 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount22;

		// dMinDist is the value used for tree-building in MAFFT
			const double dMinDist = MIN(dDist1, dDist2);
			DF.SetDist(uSeq1, uSeq2, (float) dMinDist);

			//const double dEstimatedPctId = TupleDistToEstimatedPctId(dMinDist);
			//g_dfPwId.SetDist(uSeq1, uSeq2, dEstimatedPctId);
		// **** TODO **** why does this make score slightly worse??
			//const double dKimuraDist = KimuraDist(dEstimatedPctId);
			//DF.SetDist(uSeq1, uSeq2, dKimuraDist);
			}
		}
	ProgressStepsDone();

	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		delete[] uCommonTupleCount[n];
		delete[] Letters[n];
		}
	delete[] uCommonTupleCount;
	delete[] Letters;
}
*/

inline int KmersPairwiseAlign::getElem1(int i)
{
    return (*_ptrToSeq1)[i+1];
}

inline int KmersPairwiseAlign::getElem2(int i)
{
    return (*_ptrToSeq2)[i+1];
}

}

/*******************************************************************************************************************************************************************

#include "muscle.h"
#include "distfunc.h"
#include "seqvect.h"
#include <math.h>

#define TRACE 0

#define MIN(x, y)	(((x) < (y)) ? (x) : (y))
#define MAX(x, y)	(((x) > (y)) ? (x) : (y))

const unsigned TUPLE_COUNT = 6*6*6*6*6*6;
static unsigned char Count1[TUPLE_COUNT];
static unsigned char Count2[TUPLE_COUNT];

// Nucleotide groups according to MAFFT (sextet5)
// 0 =  A
// 1 =  C
// 2 =  G
// 3 =  T
// 4 =  other

static unsigned ResidueGroup[] =
	{
	0,		// NX_A,
	1,		// NX_C,
	2,		// NX_G,
	3,		// NX_T/U
	4,		// NX_N,
	4,		// NX_R,
	4,		// NX_Y,
	4,		// NX_GAP
	};
static unsigned uResidueGroupCount = sizeof(ResidueGroup)/sizeof(ResidueGroup[0]);

static char *TupleToStr(int t)
	{
	static char s[7];
	int t1, t2, t3, t4, t5, t6;

	t1 = t%6;
	t2 = (t/6)%6;
	t3 = (t/(6*6))%6;
	t4 = (t/(6*6*6))%6;
	t5 = (t/(6*6*6*6))%6;
	t6 = (t/(6*6*6*6*6))%6;

	s[5] = '0' + t1;
	s[4] = '0' + t2;
	s[3] = '0' + t3;
	s[2] = '0' + t4;
	s[1] = '0' + t5;
	s[0] = '0' + t6;
	return s;
	}

static void ListCount(const unsigned char Count[])
	{
	for (unsigned n = 0; n < TUPLE_COUNT; ++n)
		{
		if (0 == Count[n])
			continue;
		Log("%s  %u\n", TupleToStr(n), Count[n]);
		}
	}

void DistKmer4_6(const SeqVect &v, DistFunc &DF)
	{
	if (ALPHA_DNA != g_Alpha && ALPHA_RNA != g_Alpha)
		Quit("DistKmer4_6 requires nucleo alphabet");

	const unsigned uSeqCount = v.Length();

	DF.SetCount(uSeqCount);
	if (0 == uSeqCount)
		return;

// Initialize distance matrix to zero
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			DF.SetDist(uSeq1, uSeq2, 0);
		}

// Convert to letters
	unsigned **Letters = new unsigned *[uSeqCount];
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = *(v[uSeqIndex]);
		const unsigned uSeqLength = s.Length();
		unsigned *L = new unsigned[uSeqLength];
		Letters[uSeqIndex] = L;
		for (unsigned n = 0; n < uSeqLength; ++n)
			{
			char c = s[n];
			L[n] = CharToLetterEx(c);
			if (L[n] >= 4)
				L[n] = 4;
			}
		}

	unsigned **uCommonTupleCount = new unsigned *[uSeqCount];
	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		uCommonTupleCount[n] = new unsigned[uSeqCount];
		memset(uCommonTupleCount[n], 0, uSeqCount*sizeof(unsigned));
		}

	const unsigned uPairCount = (uSeqCount*(uSeqCount + 1))/2;
	unsigned uCount = 0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &seq1 = *(v[uSeq1]);
		const unsigned uSeqLength1 = seq1.Length();
		if (uSeqLength1 < 5)
			continue;

		const unsigned uTupleCount = uSeqLength1 - 5;
		const unsigned *L = Letters[uSeq1];
		CountTuples(L, uTupleCount, Count1);
#if	TRACE
		{
		Log("Seq1=%d\n", uSeq1);
		Log("Groups:\n");
		for (unsigned n = 0; n < uSeqLength1; ++n)
			Log("%u", ResidueGroup[L[n]]);
		Log("\n");

		Log("Tuples:\n");
		ListCount(Count1);
		}
#endif

		SetProgressDesc("K-mer dist pass 1");
		for (unsigned uSeq2 = 0; uSeq2 <= uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;
			Seq &seq2 = *(v[uSeq2]);
			const unsigned uSeqLength2 = seq2.Length();
			if (uSeqLength2 < 5)
				{
				if (uSeq1 == uSeq2)
					DF.SetDist(uSeq1, uSeq2, 0);
				else
					DF.SetDist(uSeq1, uSeq2, 1);
				continue;
				}

		// First pass through seq 2 to count tuples
			const unsigned uTupleCount = uSeqLength2 - 5;
			const unsigned *L = Letters[uSeq2];
			CountTuples(L, uTupleCount, Count2);
#if	TRACE
			Log("Seq2=%d Counts=\n", uSeq2);
			ListCount(Count2);
#endif

		// Second pass to accumulate sum of shared tuples
		// MAFFT defines this as the sum over unique tuples
		// in seq2 of the minimum of the number of tuples found
		// in the two sequences.
			unsigned uSum = 0;
			for (unsigned n = 0; n < uTupleCount; ++n)
				{
				const unsigned uTuple = GetTuple(L, n);
				uSum += MIN(Count1[uTuple], Count2[uTuple]);

			// This is a hack to make sure each unique tuple counted only once.
				Count2[uTuple] = 0;
				}
#if	TRACE
			{
			Seq &s1 = *(v[uSeq1]);
			Seq &s2 = *(v[uSeq2]);
			const char *pName1 = s1.GetName();
			const char *pName2 = s2.GetName();
			Log("Common count %s(%d) - %s(%d) =%u\n",
			  pName1, uSeq1, pName2, uSeq2, uSum);
			}
#endif
			uCommonTupleCount[uSeq1][uSeq2] = uSum;
			uCommonTupleCount[uSeq2][uSeq1] = uSum;
			}
		}
	ProgressStepsDone();

	uCount = 0;
	SetProgressDesc("K-mer dist pass 2");
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &s1 = *(v[uSeq1]);
		const char *pName1 = s1.GetName();

		double dCommonTupleCount11 = uCommonTupleCount[uSeq1][uSeq1];
		if (0 == dCommonTupleCount11)
			dCommonTupleCount11 = 1;

		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;

			double dCommonTupleCount22 = uCommonTupleCount[uSeq2][uSeq2];
			if (0 == dCommonTupleCount22)
				dCommonTupleCount22 = 1;

			const double dDist1 = 3.0*(dCommonTupleCount11 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount11;
			const double dDist2 = 3.0*(dCommonTupleCount22 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount22;

		// dMinDist is the value used for tree-building in MAFFT
			const double dMinDist = MIN(dDist1, dDist2);
			DF.SetDist(uSeq1, uSeq2, (float) dMinDist);

			//const double dEstimatedPctId = TupleDistToEstimatedPctId(dMinDist);
			//g_dfPwId.SetDist(uSeq1, uSeq2, dEstimatedPctId);
		// **** TODO **** why does this make score slightly worse??
			//const double dKimuraDist = KimuraDist(dEstimatedPctId);
			//DF.SetDist(uSeq1, uSeq2, dKimuraDist);
			}
		}
	ProgressStepsDone();

	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		delete[] uCommonTupleCount[n];
		delete[] Letters[n];
		}
	delete[] uCommonTupleCount;
	delete[] Letters;
	}


#include "muscle.h"
#include "distfunc.h"
#include "seqvect.h"
#include <math.h>

#define TRACE 0

#define MIN(x, y)	(((x) < (y)) ? (x) : (y))
#define MAX(x, y)	(((x) > (y)) ? (x) : (y))

const unsigned TUPLE_COUNT = 6*6*6*6*6*6;
static unsigned char Count1[TUPLE_COUNT];
static unsigned char Count2[TUPLE_COUNT];

// Amino acid groups according to MAFFT (sextet5)
// 0 =  A G P S T
// 1 =  I L M V
// 2 =  N D Q E B Z
// 3 =  R H K
// 4 =  F W Y
// 5 =  C
// 6 =  X . - U
unsigned ResidueGroup[] =
	{
	0,		// AX_A,
	5,		// AX_C,
	2,		// AX_D,
	2,		// AX_E,
	4,		// AX_F,
	0,		// AX_G,
	3,		// AX_H,
	1,		// AX_I,
	3,		// AX_K,
	1,		// AX_L,
	1,		// AX_M,
	2,		// AX_N,
	0,		// AX_P,
	2,		// AX_Q,
	3,		// AX_R,
	0,		// AX_S,
	0,		// AX_T,
	1,		// AX_V,
	4,		// AX_W,
	4,		// AX_Y,

	2,		// AX_B,	// D or N
	2,		// AX_Z,	// E or Q
	0,		// AX_X,	// Unknown		// ******** TODO *************
										// This isn't the correct way of avoiding group 6
	0		// AX_GAP,					// ******** TODO ******************
	};
unsigned uResidueGroupCount = sizeof(ResidueGroup)/sizeof(ResidueGroup[0]);

static char *TupleToStr(int t)
	{
	static char s[7];
	int t1, t2, t3, t4, t5, t6;

	t1 = t%6;
	t2 = (t/6)%6;
	t3 = (t/(6*6))%6;
	t4 = (t/(6*6*6))%6;
	t5 = (t/(6*6*6*6))%6;
	t6 = (t/(6*6*6*6*6))%6;

	s[5] = '0' + t1;
	s[4] = '0' + t2;
	s[3] = '0' + t3;
	s[2] = '0' + t4;
	s[1] = '0' + t5;
	s[0] = '0' + t6;
	return s;
	}

static unsigned GetTuple(const unsigned uLetters[], unsigned n)
	{
	assert(uLetters[n] < uResidueGroupCount);
	assert(uLetters[n+1] < uResidueGroupCount);
	assert(uLetters[n+2] < uResidueGroupCount);
	assert(uLetters[n+3] < uResidueGroupCount);
	assert(uLetters[n+4] < uResidueGroupCount);
	assert(uLetters[n+5] < uResidueGroupCount);

	unsigned u1 = ResidueGroup[uLetters[n]];
	unsigned u2 = ResidueGroup[uLetters[n+1]];
	unsigned u3 = ResidueGroup[uLetters[n+2]];
	unsigned u4 = ResidueGroup[uLetters[n+3]];
	unsigned u5 = ResidueGroup[uLetters[n+4]];
	unsigned u6 = ResidueGroup[uLetters[n+5]];

	return u6 + u5*6 + u4*6*6 + u3*6*6*6 + u2*6*6*6*6 + u1*6*6*6*6*6;
	}

static void CountTuples(const unsigned L[], unsigned uTupleCount, unsigned char Count[])
	{
	memset(Count, 0, TUPLE_COUNT*sizeof(unsigned char));
	for (unsigned n = 0; n < uTupleCount; ++n)
		{
		const unsigned uTuple = GetTuple(L, n);
		++(Count[uTuple]);
		}
	}

static void ListCount(const unsigned char Count[])
	{
	for (unsigned n = 0; n < TUPLE_COUNT; ++n)
		{
		if (0 == Count[n])
			continue;
		Log("%s  %u\n", TupleToStr(n), Count[n]);
		}
	}

void DistKmer6_6(const SeqVect &v, DistFunc &DF)
	{
	const unsigned uSeqCount = v.Length();

	DF.SetCount(uSeqCount);
	if (0 == uSeqCount)
		return;

// Initialize distance matrix to zero
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			DF.SetDist(uSeq1, uSeq2, 0);
		}

// Convert to letters
	unsigned **Letters = new unsigned *[uSeqCount];
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = *(v[uSeqIndex]);
		const unsigned uSeqLength = s.Length();
		unsigned *L = new unsigned[uSeqLength];
		Letters[uSeqIndex] = L;
		for (unsigned n = 0; n < uSeqLength; ++n)
			{
			char c = s[n];
			L[n] = CharToLetterEx(c);
			assert(L[n] < uResidueGroupCount);
			}
		}

	unsigned **uCommonTupleCount = new unsigned *[uSeqCount];
	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		uCommonTupleCount[n] = new unsigned[uSeqCount];
		memset(uCommonTupleCount[n], 0, uSeqCount*sizeof(unsigned));
		}

	const unsigned uPairCount = (uSeqCount*(uSeqCount + 1))/2;
	unsigned uCount = 0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &seq1 = *(v[uSeq1]);
		const unsigned uSeqLength1 = seq1.Length();
		if (uSeqLength1 < 5)
			continue;

		const unsigned uTupleCount = uSeqLength1 - 5;
		const unsigned *L = Letters[uSeq1];
		CountTuples(L, uTupleCount, Count1);
#if	TRACE
		{
		Log("Seq1=%d\n", uSeq1);
		Log("Groups:\n");
		for (unsigned n = 0; n < uSeqLength1; ++n)
			Log("%u", ResidueGroup[L[n]]);
		Log("\n");

		Log("Tuples:\n");
		ListCount(Count1);
		}
#endif

		SetProgressDesc("K-mer dist pass 1");
		for (unsigned uSeq2 = 0; uSeq2 <= uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;
			Seq &seq2 = *(v[uSeq2]);
			const unsigned uSeqLength2 = seq2.Length();
			if (uSeqLength2 < 5)
				{
				if (uSeq1 == uSeq2)
					DF.SetDist(uSeq1, uSeq2, 0);
				else
					DF.SetDist(uSeq1, uSeq2, 1);
				continue;
				}

		// First pass through seq 2 to count tuples
			const unsigned uTupleCount = uSeqLength2 - 5;
			const unsigned *L = Letters[uSeq2];
			CountTuples(L, uTupleCount, Count2);
#if	TRACE
			Log("Seq2=%d Counts=\n", uSeq2);
			ListCount(Count2);
#endif

		// Second pass to accumulate sum of shared tuples
		// MAFFT defines this as the sum over unique tuples
		// in seq2 of the minimum of the number of tuples found
		// in the two sequences.
			unsigned uSum = 0;
			for (unsigned n = 0; n < uTupleCount; ++n)
				{
				const unsigned uTuple = GetTuple(L, n);
				uSum += MIN(Count1[uTuple], Count2[uTuple]);

			// This is a hack to make sure each unique tuple counted only once.
				Count2[uTuple] = 0;
				}
#if	TRACE
			{
			Seq &s1 = *(v[uSeq1]);
			Seq &s2 = *(v[uSeq2]);
			const char *pName1 = s1.GetName();
			const char *pName2 = s2.GetName();
			Log("Common count %s(%d) - %s(%d) =%u\n",
			  pName1, uSeq1, pName2, uSeq2, uSum);
			}
#endif
			uCommonTupleCount[uSeq1][uSeq2] = uSum;
			uCommonTupleCount[uSeq2][uSeq1] = uSum;
			}
		}
	ProgressStepsDone();

	uCount = 0;
	SetProgressDesc("K-mer dist pass 2");
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &s1 = *(v[uSeq1]);
		const char *pName1 = s1.GetName();

		double dCommonTupleCount11 = uCommonTupleCount[uSeq1][uSeq1];
		if (0 == dCommonTupleCount11)
			dCommonTupleCount11 = 1;

		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;

			double dCommonTupleCount22 = uCommonTupleCount[uSeq2][uSeq2];
			if (0 == dCommonTupleCount22)
				dCommonTupleCount22 = 1;

			const double dDist1 = 3.0*(dCommonTupleCount11 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount11;
			const double dDist2 = 3.0*(dCommonTupleCount22 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount22;

		// dMinDist is the value used for tree-building in MAFFT
			const double dMinDist = MIN(dDist1, dDist2);
			DF.SetDist(uSeq1, uSeq2, (float) dMinDist);

			//const double dEstimatedPctId = TupleDistToEstimatedPctId(dMinDist);
			//g_dfPwId.SetDist(uSeq1, uSeq2, dEstimatedPctId);
		// **** TODO **** why does this make score slightly worse??
			//const double dKimuraDist = KimuraDist(dEstimatedPctId);
			//DF.SetDist(uSeq1, uSeq2, dKimuraDist);
			}
		}
	ProgressStepsDone();

	for (unsigned n = 0; n < uSeqCount; ++n)
		delete[] uCommonTupleCount[n];
	delete[] uCommonTupleCount;
	delete[] Letters;
	}

double PctIdToMAFFTDist(double dPctId)
	{
	if (dPctId < 0.05)
		dPctId = 0.05;
	double dDist = -log(dPctId);
	return dDist;
	}

double PctIdToHeightMAFFT(double dPctId)
	{
	return PctIdToMAFFTDist(dPctId);
	}
*/
