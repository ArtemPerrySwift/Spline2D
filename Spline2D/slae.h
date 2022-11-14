#pragma once
#pragma once
#include "sparsematrix.h"

using namespace matrix;
typedef double real;
namespace slae
{
	template<class MatrixType, typename = std::enable_if_t<std::is_base_of<SparseMatrix, MatrixType>::value>>
	struct SLAE
	{
		MatrixType A;
		double* b;
		void init(MatrixType A, double* b)
		{
			this->A = A;
			this->b = b;
		}

		void init(MatrixType A)
		{
			b = new double[A.n];
			init(A, b);
		}

		void setOneVariableSolve(int iVar, double varMean);
	};
	int solveSLAU3(double A[3][3], double b[3], double x[3]);

}