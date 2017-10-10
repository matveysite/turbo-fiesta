#pragma once
#include "Matrix.h"
#include "Set.h"

void dualSimplexMethod(Matrix& A, Matrix& b, Matrix& c, Matrix& dl, Matrix& dh, Set& base);
void simplexMethod(Matrix& A, Matrix& b, Matrix& c, Matrix& dl, Matrix& dh, Matrix& x, Set& base);
void problemSensitivity(Matrix& A, Matrix& b, Matrix& c, Matrix& dl, Matrix& dh, Set& optimalBase);
double sign(double);