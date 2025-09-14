/**
 * Author: Mario Joao Jr.
 * 
 */

#ifndef CONSTRENTRY_H
#define CONSTRENTRY_H

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>


using namespace std;

namespace clustalw
{

class ConstrEntry
{
    public:
        /* Functions */
	ConstrEntry(int s1, int s2, int r1, int r2, int weight, int cons, int ind);

	int getS1() {return s1;};
	int getS2() {return s2;};
	int getR1() {return r1;};
	int getR2() {return r2;};
	int getWeight() {return weight;};
	int getCons() {return cons;};
	int getInd() {return ind;};
	void setInd(int pind) {ind = pind;};
        
/*
        const vector<int>* getSequence(int index){return &seqArray[index];}; // For Pairwise!
        const vector<int>* getSequence(int index) const {return &seqArray[index];};
        const vector<int>* getSequenceFromUniqueId(unsigned long id); // For iteration        
        const SeqArray* getSeqArray() const {return &seqArray;}; // For multiple align!
        SeqArray* getSeqArrayForRealloc(){return &seqArray;};
*/
        /* Attributes */
	// s1 always < s2 
	int	s1;
	int	s2;
	int	r1;
	int	r2;
	int	weight;
	int	cons;
	int	ind;

    private:
        /* Functions */
        
        /* Attributes */
};
}
#endif

