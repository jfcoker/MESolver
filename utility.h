#pragma once
#include "pch.h"
#include "site.h"

void printMatrix(gsl_matrix* m);
void printVector(gsl_vector* v, bool horizontal = false);
void printDiagonal(gsl_matrix* m, bool horizontal = false);
void printOccProbs(std::vector<site>& sites, int precision = 2);