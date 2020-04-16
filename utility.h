#pragma once
#include "pch.h"

void printMatrix(gsl_matrix* m);
void printVector(gsl_vector* v, bool horizontal = true);
void printDiagonal(gsl_matrix* m, bool horizontal = true);