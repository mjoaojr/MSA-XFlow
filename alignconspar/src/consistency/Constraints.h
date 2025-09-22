/**
 * Author: Mario Joao Jr.
 * 
 */

#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#ifndef CONSTRENTRY_H
	#include "ConstrEntry.h"
#endif

#ifndef EDGE_H
	#include "Edge.h"
#endif

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <unordered_map>

//#include "SparseMatrix.h"

#define	KEY(x,y)	(((unsigned long)(x)<<32)+(y))

using namespace std;

struct ExtEdge
{
/*
    public:
         Functions 
		ExtEdge (void) {p1 = 0; p2=0; diag=0; score=0;}
		ExtEdge (int a, int b, int c, float d) {p1 = a; p2=b; diag=c; score=d;}
*/
        /* Attributes */
	int	p1,
		p2;
//		diag;
	float	score;
};

namespace clustalw
{

class Constraints
{
    public:
        /* Functions */
        Constraints(int nSeqs);
	void addEntry(ConstrEntry entry);
	int getIniCLPair (int s1, int s2);
	void print (void);
	void printEdges (void);
	void printResIndex (void);
	void saveEdges (string filename);
	void loadEdges (string fileName);

/*
	void setSparseMatrix (int a, int b, int lenA, int lenB, vector<float> post){sparseMatrices[a][b] = new SparseMatrix (lenA, lenB, post); sparseMatrices[b][a] = NULL;}
	vector<vector<SparseMatrix *>>getSparseMatrices (){return sparseMatrices;}
	void setSparseMatrices (vector<vector<SparseMatrix *>> sm){sparseMatrices = sm;}
*/
	vector<vector<vector<Edge>>> *getEdges(){return &resEdges;}
	int getMaxExtValue () {return maxExtValue;}
	void setMaxExtValue (int v) {maxExtValue = v;}
	int getMaxValue () {return maxValue;}
	void setMaxValue (int v) {maxValue = v;}
	int getNEdges () {return nEdges;}
	void setNEdges (int e) {nEdges = e;}
	int getNoMatch () {return noMatch;}
	void setNoMatch (int m) {noMatch = m;}
	int getNormalise () {return normalise;}
	void setNormalise (int m) {normalise = m;}
        
        /* Attributes */
	int nSeqs;
	vector<ConstrEntry> constrList;
	vector<vector<int>> resIndex;
	vector<vector<vector<Edge>>> resEdges;
	int nEdges;
//	vector<vector<SparseMatrix *>>sparseMatrices;
	int maxValue;
	int maxExtValue;
	int noMatch;
	int normalise;
	vector<ExtEdge> *extEdges;
	vector<int> 	*indExtEdges;
	std::unordered_map<unsigned long int, float>	*umap;

    private:
        /* Functions */
	void moveRight (int ini);
	void addEntry2 (ConstrEntry entry);
	void addEdge (ConstrEntry entry);
	void setEdge (vector<Edge> &vet, int s2, int r2, int weight);
        
        void clearSeqArray();
        /* Attributes */
};
}
#endif

