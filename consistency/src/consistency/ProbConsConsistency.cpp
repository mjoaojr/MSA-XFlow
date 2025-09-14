/**
 * Author: Mario Joao Jr.
 * 
 * Adapted from ProbCons Code
 */

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>

#include "ProbConsConsistency.h"


/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the consistency transformation.  The
// formula used is:
//                     1
//    P'(x[i]-y[j]) = ---  sum   sum P(x[i]-z[k]) P(z[k]-y[j])
//                    |S| z in S  k
//
// where S = {x, y, all other sequences...}
//
/////////////////////////////////////////////////////////////////

vector<vector<SparseMatrix *> > DoRelaxation (clustalw::Alignment *alignPtr, vector<vector<SparseMatrix *> > &sparseMatrices)
{
	const int numSeqs = alignPtr->getNumSeqs();
	vector<vector<SparseMatrix *> > newSparseMatrices (numSeqs, vector<SparseMatrix *>(numSeqs, NULL));

	cerr << endl;
	// for every pair of sequences
	for (int i = 0; i < numSeqs; i++)
	{
		for (int j = i+1; j < numSeqs; j++)
		{
			const vector<int> *seq1 = alignPtr->getSequence (i+1);
			const vector<int> *seq2 = alignPtr->getSequence (j+1);

//			if (enableVerbose)
//			cerr << "Relaxing (" << i+1 << ") " << seq1->GetHeader() << " vs. " << "(" << j+1 << ") " << seq2->GetHeader() << ": ";
//			cerr << "Relaxing (" << i+1 << ") " << " vs. " << "(" << j+1 << ") " << ": ";

			// get the original posterior matrix
			vector<float> *posteriorPtr = sparseMatrices[i][j]->GetPosterior();
			vector<float> &posterior = *posteriorPtr;

			const int seq1Length = alignPtr->getSeqLength(i + 1);
			const int seq2Length = alignPtr->getSeqLength(j + 1);

			// contribution from the summation where z = x and z = y
			for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] += posterior[k];

//			if (enableVerbose)
//			cerr << sparseMatrices[i][j]->GetNumCells() << " --> ";

			// contribution from all other sequences
			for (int k = 0; k < numSeqs; k++) 
			{
				if (k != i && k != j)
				{
					if (k < i)
						Relax1 (sparseMatrices[k][i], sparseMatrices[k][j], posterior);
					else 
						if (k > i && k < j)
							Relax (sparseMatrices[i][k], sparseMatrices[k][j], posterior);
						else 
						{
							SparseMatrix *temp = sparseMatrices[j][k]->ComputeTranspose();
							Relax (sparseMatrices[i][k], temp, posterior);
							delete temp;
						}
				}
			}


			// now renormalization
			for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) 
				posterior[k] /= numSeqs;

			// mask out positions not originally in the posterior matrix
			SparseMatrix *matXY = sparseMatrices[i][j];
			for (int y = 0; y <= seq2Length; y++) 
				posterior[y] = 0;
			for (int x = 1; x <= seq1Length; x++)
			{
				vector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
				vector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
				vector<float>::iterator base = posterior.begin() + x * (seq2Length + 1);
				int curr = 0;
				while (XYptr != XYend)
				{
					// zero out all cells until the first filled column
					while (curr < XYptr->first)
					{
						base[curr] = 0;
						curr++;
					}

					// now, skip over this column
					curr++;
					++XYptr;
				}

				// zero out cells after last column
				while (curr <= seq2Length)
				{
					base[curr] = 0;
					curr++;
				}
			}

			// save the new posterior matrix
			newSparseMatrices[i][j] = new SparseMatrix (seq1Length, seq2Length, posterior);
			newSparseMatrices[j][i] = NULL;

//			if (enableVerbose)
//			cerr << newSparseMatrices[i][j]->GetNumCells() << " -- ";

			delete posteriorPtr;

//			if (enableVerbose)
//			cerr << "done." << endl;
		}
	}

	return newSparseMatrices;
}


/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, vector<float> &posterior)
{
	int lengthX = matXZ->GetSeq1Length();
	int lengthY = matZY->GetSeq2Length();

	// for every x[i]
	for (int i = 1; i <= lengthX; i++)
	{
		vector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
		vector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

		vector<float>::iterator base = posterior.begin() + i * (lengthY + 1);

		// iterate through all x[i]-z[k]
		while (XZptr != XZend)
		{
			vector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
			vector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
			const float XZval = XZptr->second;

			// iterate through all z[k]-y[j]
			while (ZYptr != ZYend)
			{
				base[ZYptr->first] += XZval * ZYptr->second;
				ZYptr++;
			}
			XZptr++;
		}
	}
}

/////////////////////////////////////////////////////////////////
// Relax1()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax1 (SparseMatrix *matZX, SparseMatrix *matZY, vector<float> &posterior)
{
	int lengthZ = matZX->GetSeq1Length();
	int lengthY = matZY->GetSeq2Length();

	// for every z[k]
	for (int k = 1; k <= lengthZ; k++)
	{
		vector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
		vector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

		// iterate through all z[k]-x[i]
		while (ZXptr != ZXend)
		{
			vector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
			vector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
			const float ZXval = ZXptr->second;
			vector<float>::iterator base = posterior.begin() + ZXptr->first * (lengthY + 1);

			// iterate through all z[k]-y[j]
			while (ZYptr != ZYend)
			{
				base[ZYptr->first] += ZXval * ZYptr->second;
				ZYptr++;
			}
			ZXptr++;
		}
	}
}

