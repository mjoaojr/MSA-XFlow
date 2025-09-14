/**
 * Author: Mario Joao Jr.
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "LcsPairwiseAlign.h"
#include <math.h>

namespace clustalw
{

LcsPairwiseAlign::LcsPairwiseAlign()
{
	_maxAlnLength = 0;
}

void LcsPairwiseAlign::pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, int iEnd, int jStart, int jEnd)
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
    
        for (si = utilityObject->MAX(0, iStart); si < _numSeqs && si < iEnd; si++)
        {
            n = alignPtr->getSeqLength(si + 1);

            for (sj = utilityObject->MAX(si+1, jStart+1); sj < _numSeqs && sj < jEnd; sj++)
            {
                m = alignPtr->getSeqLength(sj + 1);
                // align the sequences
            
                _ptrToSeq1 = alignPtr->getSequence(si + 1);
                _ptrToSeq2 = alignPtr->getSequence(sj + 1);
            
		// Use LCS to calculate score
		lcsScore = algc (0, n-1, 0, m-1);
//		lcsScore = min(m, n) - lcsScore;
		_score = (100.0 * (float)lcsScore) / (float)(min(m, n));

                _score = ((float)100.0 - _score) / (float)100.0;
                distMat->SetAt(si + 1, sj + 1, _score);
                distMat->SetAt(sj + 1, si + 1, _score);
                
                if(userParameters->getDisplayInfo())
                {
                    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %.4f",
                                        si+1, sj+1, (float)_score);     
                }
            }
        }
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the LcsPairwiseAlign class.\n"
             << e.what() << "\n";
        exit(1);    
    }
}

void LcsPairwiseAlign::algb (int inia, int fima, int inib, int fimb, vector<TELEM> &LL)
{
	int	i,
		j,
		m,
		n;

	m = fima - inia + 1;
	n = fimb - inib + 1;

	vector<TELEM> K0(n+1, 0);
	vector<TELEM> K1(n+1);
        K1[0] = 0;

	for (i = inia; i <= fima; i++)
	{
                if ((i - inia) % 2 == 0)
                {
			for (j = inib; j <= fimb; j++)
			{
				if (getElem1(i) == getElem2(j)) 
					K1[j-inib+1] = 1 + K0[j-inib];
				else 
					K1[j-inib+1] = max(K0[j-inib+1], K1[j-inib]);
			}
		}
		else
		{
			for (j = inib; j <= fimb; j++)
			{
				if (getElem1(i) == getElem2(j)) 
					K0[j-inib+1] = 1 + K1[j-inib];
				else 
					K0[j-inib+1] = max(K1[j-inib+1], K0[j-inib]);
			}
		}
	}

        if (m % 2)
	{
//		memcpy (LL, K1+1, n * sizeof (TELEM));
		K1.erase(K1.begin());
		LL = K1;
	}
	else
	{
//		memcpy (LL, K0+1, n * sizeof (TELEM));
		K0.erase(K0.begin());
		LL = K0;
	}

	K0.clear();
	K1.clear();
}

void LcsPairwiseAlign::algbback (int inia, int fima, int inib, int fimb, vector<TELEM> &LL)
{
	int	i,
		j,
		m,
		n;

	m = fima - inia + 1;
	n = fimb - inib + 1;

	vector<TELEM> K0(n+1, 0);
	vector<TELEM> K1(n+1);
        K1[n] = 0;

	for (i = fima; i >= inia; i--)
	{
                if ((fima - i) % 2 == 0)
                {
			for (j = fimb; j >= inib; j--)
			{
				if (getElem1(i) == getElem2(j)) 
					K1[j-inib] = 1 + K0[j-inib+1];
				else 
					K1[j-inib] = max(K0[j-inib], K1[j-inib+1]);
			}
		}
		else
		{
			for (j = fimb; j >= inib; j--)
			{
				if (getElem1(i) == getElem2(j)) 
					K0[j-inib] = 1 + K1[j-inib+1];
				else 
					K0[j-inib] = max(K1[j-inib], K0[j-inib+1]);
			}
		}
	}

        if (m % 2)
	{
//		memcpy (LL, K1, n * sizeof (TELEM));
//		K1.erase(K1.end());
		K1.pop_back();
		LL = K1;
	}
	else
	{
//		memcpy (LL, K0, n * sizeof (TELEM));
//		K0.erase(K0.end());
		K0.pop_back();
		LL = K0;
	}
		
	K0.clear();
	K1.clear();
}

int LcsPairwiseAlign::algc (int inia, int fima, int inib, int fimb)
{
	int	i,
		j,
		k,
		kk,
		m,
		meio,
		n,
		ret;
	int	C1,
		C2;
	char	c;

	m = fima - inia + 1;
	n = fimb - inib + 1;

/*
printf ("inia = %d, fima = %d, m=%d\n", inia, fima, m);
printf ("inib = %d, fimb = %d, n=%d\n", inib, fimb, n);
*/
	if ((m <= 0) || (n <= 0))
	{
//		C[0] = 0;
		return 0;
	}

	if (m == 1)
	{
//		strcpy (C, "");
		for (j = inib; j <= fimb; j++)
		{
			if (getElem1(inia) == getElem2(j))
			{
//				C[0] = a[inia];
//				C[1] = '\0';
//				break;
				return 1;
			}
		}

		return 0;
	}

	meio = inia + m/2;

	vector<TELEM> L1(n);

//printf ("Chamando1 algb com (%d,%d) e (%d,%d)\n", inia, meio-1, inib, fimb);
	algb (inia, meio-1, inib, fimb, L1);
	//prints (L1, n);

	vector<TELEM> L2(n);

//printf ("Chamando2 algb com (%d,%d) e (%d,%d)\n", meio, fima, inib, fimb);
	algbback (meio, fima, inib, fimb, L2);
	//prints (L2, n);

	for (i = 1, k = 0, kk = L2[0]; i <= fimb - inib; i++)
	{
		if (L1[i-1] + L2[i] > kk)
		{
			k = i;
			kk = L1[i-1] + L2[i];
		}
	}

	if (L1[i-1] > kk)
	{
		k = i;
		kk = L1[i];
	}

	L1.clear();
	L2.clear();

//printf ("meio = %d, k = %d\n", meio, k);
//printf ("Chamando1 algc com (%d,%d) e (%d,%d)\n", inia, meio-1, inib, inib+k-1);
	C1 = algc (inia, meio-1, inib, inib+k-1);
//printf ("Chamando2 algc com (%d,%d) e (%d,%d)\n", meio, fima, inib+k, fimb);
	C2 = algc (meio, fima, inib+k, fimb);

	return C1 + C2;
}

inline int LcsPairwiseAlign::getElem1(int i)
{
    return (*_ptrToSeq1)[i+1];
}

inline int LcsPairwiseAlign::getElem2(int i)
{
    return (*_ptrToSeq2)[i+1];
}

}
