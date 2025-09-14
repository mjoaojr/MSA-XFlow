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
#include "ConstrEntry.h"

using namespace std;

namespace clustalw
{

ConstrEntry::ConstrEntry(int s1, int s2, int r1, int r2, int weight, int cons, int ind)
{
	if (s1 < s2)
	{
		this->s1 = s1;
		this->s2 = s2;
		this->r1 = r1;
		this->r2 = r2;
	}
	else
	{
		this->s1 = s2;
		this->s2 = s1;
		this->r1 = r2;
		this->r2 = r1;
	}

	this->weight = weight;
	this->cons = cons;
	this->ind = ind;
}

}

