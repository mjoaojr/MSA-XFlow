/**
 * Author: Mario Joao Jr.
 * 
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <exception>
#include <cmath>
#include <sstream>
#include "Edge.h"

using namespace std;

namespace clustalw
{

Edge::Edge(int s2, int r2, int weight)
{
	this->s2 = s2;
	this->r2 = r2;
	this->weight = weight;
}

Edge::Edge(void)
{
	this->s2 = 0;
	this->r2 = 0;
	this->weight = 0;
}

}

