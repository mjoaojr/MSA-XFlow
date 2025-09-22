/**
 * Author: Mario Joao Jr.
 * 
 */

#ifndef EDGE_H
#define EDGE_H

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>


using namespace std;

namespace clustalw
{

class Edge
{
    public:
        /* Functions */
	Edge (int s2, int r2, int weight);
	Edge (void);

	int getS2() {return s2;};
	int getR2() {return r2;};
	int getWeight() {return weight;};
	void setEdge(int ps2, int pr2, int pweight) {s2=ps2; r2=pr2; weight=pweight;};
	void setWeight(int pweight) {weight=pweight;};
        
        /* Attributes */
	int	s2;
	int	r2;
	int	weight;

    private:
        /* Functions */
        
        /* Attributes */
};
}
#endif

