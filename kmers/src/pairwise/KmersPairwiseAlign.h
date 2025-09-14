/**
 * Author: Mario Joao Jr.
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 */
#ifndef KMERSPAIRWISEALIGN_H
#define KMERSPAIRWISEALIGN_H

#include "PairwiseAlignBase.h"

#define MAX_CHAR 256

/***
From Bioperl docs:
Extended DNA / RNA alphabet
------------------------------------------
Symbol       Meaning      Nucleic Acid
------------------------------------------
    A            A           Adenine
    C            C           Cytosine
    G            G           Guanine
    T            T           Thymine
    U            U           Uracil
    M          A or C
    R          A or G
    W          A or T
    S          C or G
    Y          C or T
    K          G or T
    V        A or C or G
    H        A or C or T
    D        A or G or T
    B        C or G or T
    X      G or A or T or C
    N      G or A or T or C

IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
         Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.
***/

// NX=Nucleotide alphabet with extensions
enum NX
{
	NX_A,
	NX_C,
	NX_G,
	NX_T,
	NX_U = NX_T,

	NX_M, // AC
	NX_R, // AG
	NX_W, // AT
	NX_S, // CG
	NX_Y, // CT
	NX_K, // GT
	NX_V, // ACG
	NX_H, // ACT
	NX_D, // AGT
	NX_B, // CGT
	NX_X, // GATC
	NX_N, // GATC
	NX_GAP
};

// AX=Amino alphabet with eXtensions (B, Z and X)
enum AX
{
        AX_A,
        AX_B,   // D or N
        AX_C,
        AX_D,
        AX_E,
        AX_F,
        AX_G,
        AX_H,
        AX_I,
        AX_K,
        AX_L,
        AX_M,
        AX_N,
        AX_O,   // Any
        AX_P,
        AX_Q,
        AX_R,
        AX_S,
        AX_T,
        AX_U,   // Any
        AX_V,
        AX_W,
        AX_X,   // Any
        AX_Y,
        AX_Z,   // E or Q

        AX_GAP,
};



namespace clustalw
{

class KmersPairwiseAlign : public PairwiseAlignBase
{
    public:
        /* Functions */
        KmersPairwiseAlign();
	virtual ~KmersPairwiseAlign(){};

        virtual void pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, 
                                   int iEnd, int jStart, int jEnd); 
        /* Attributes */

    private:
        /* Functions */
	int getElem1(int i);
	int getElem2(int i);
/*
	int algc (int inia, int fima, int inib, int fimb);
	void algb (int inia, int fima, int inib, int fimb, vector<TELEM> &LL);
	void algbback (int inia, int fima, int inib, int fimb, vector<TELEM> &LL);
*/
	void CountTuples(const unsigned L[], unsigned uTupleCount, unsigned char Count[]);
	unsigned GetTuple(const unsigned uLetters[], unsigned n);

        /* Attributes */
        // I have constant pointers to the data. This allows for the fastest access.
        const vector<int>* _ptrToSeq1;
        const vector<int>* _ptrToSeq2;
        int _maxAlnLength;
	int	kmersScore;
	vector <unsigned> CharToLetterEx = vector<unsigned>(MAX_CHAR);

	static const unsigned TUPLE_COUNT = 6*6*6*6*6*6;
	unsigned char Count1[TUPLE_COUNT];
	unsigned char Count2[TUPLE_COUNT];

	// Nucleotide groups according to MAFFT (sextet5)
	// 0 =  A
	// 1 =  C
	// 2 =  G
	// 3 =  T
	// 4 =  other

	const vector<unsigned char> DNAResidueGroup =
	{
		0,              // NX_A,
		1,              // NX_C,
		2,              // NX_G,
		3,              // NX_T/U
		4,              // NX_N,
		4,              // NX_R,
		4,              // NX_Y,
		4,              // NX_GAP
	};
	const unsigned uDNAResidueGroupCount = sizeof(ResidueGroup)/sizeof(ResidueGroup[0]);

	// Amino acid groups according to MAFFT (sextet5)
	// 0 =  A G P S T
	// 1 =  I L M V
	// 2 =  N D Q E B Z
	// 3 =  R H K
	// 4 =  F W Y
	// 5 =  C
	// 6 =  X . - U
	const vector<unsigned char> ProtResidueGroup =
	{
		0,              // AX_A,
		5,              // AX_C,
		2,              // AX_D,
		2,              // AX_E,
		4,              // AX_F,
		0,              // AX_G,
		3,              // AX_H,
		1,              // AX_I,
		3,              // AX_K,
		1,              // AX_L,
		1,              // AX_M,
		2,              // AX_N,
		0,              // AX_P,
		2,              // AX_Q,
		3,              // AX_R,
		0,              // AX_S,
		0,              // AX_T,
		1,              // AX_V,
		4,              // AX_W,
		4,              // AX_Y,

		2,              // AX_B,        // D or N
		2,              // AX_Z,        // E or Q
		0,              // AX_X,        // Unknown              // ******** TODO *************
											// This isn't the correct way of avoiding group 6
		0               // AX_GAP,                                      // ******** TODO ******************
	};
	unsigned uProtResidueGroupCount = sizeof(ResidueGroup)/sizeof(ResidueGroup[0]);


	vector<unsigned char> ResidueGroup;
	unsigned uResidueGroupCount;

/*
        int intScale;
        float mmScore;
        int printPtr;
        int lastPrint;

        int _gapOpen; // scaled to be an integer, this is not a mistake
        int _gapExtend; // scaled to be an integer, not a mistake
        int seq1;
        int seq2;
        int matrix[NUMRES][NUMRES];
        int maxScore;
        int sb1;
        int sb2;
        int se1;
        int se2;
*/

};

}
#endif
