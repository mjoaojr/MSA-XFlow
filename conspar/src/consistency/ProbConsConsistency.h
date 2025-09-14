/**
 * Author: Mario Joao Jr.
 * 
 * Adapted from ProbCons Code
 */

#ifndef PROBCONSCONS_H
#define PROBCONSCONS_H

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>

#include "SparseMatrix.h"
#include "../alignment/Alignment.h"

vector<vector<SparseMatrix *> > DoRelaxation (clustalw::Alignment *alignPtr, vector<vector<SparseMatrix *> > &sparseMatrices);
void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, vector<float> &posterior);
void Relax1 (SparseMatrix *matZX, SparseMatrix *matZY, vector<float> &posterior);

#endif

