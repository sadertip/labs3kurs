#include "Gaus.cpp"
#include "Matrices_template.cpp"
#include "QR decomposition.h"
#include "other.cpp"


int main()
{
 using T = double;
	matrix<T> A = matrix<T>(3, 3, { 1,0,4,0,1,1,1,2,1 });
	matrix<T> b = matrix<T>(3, 1, { 1, 2, 3 });
	QR_decomposition(A, b);
	return 0;
}