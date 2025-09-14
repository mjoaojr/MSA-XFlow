/**
 * Author: Mario Joao Jr.
 * 
 * Adapted from TCoffee Code
 */

#ifndef TCOFFEECONS_H
#define TCOFFEECONS_H

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>

#include "../alignment/Alignment.h"

void relaxCL (clustalw::Alignment *alignPtr);
void maxConstraints (clustalw::Alignment *alignPtr);
void cl2PairListExt (clustalw::Alignment *alignPtr, vector<int> *group);


#endif

