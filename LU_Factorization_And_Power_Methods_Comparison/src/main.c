#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#pragma warning (disable : 4996)
#define SIZE 12

unsigned int luOps1, luOps2;
unsigned int itOps1, itOps2;

double NormMatrix(double* matr, int size);
double* InvMatrix(double* matr, int size);
double GetCond(double* matr, int size);

double* CreateMatrix(int size) {
	int i, j, a = 0;
	double* matr = (double*)malloc(sizeof(double) * size * size), sum;
	
	for (i = 0; i < size; i++) {
		sum = 0;
		for (j = 0; j < size; j++) {
			a = rand() % 9 + 1;
			matr[i * size + j] = (double)(a);
			sum += (double)a;
		}
		matr[i * size + i] = sum + 5;
	}

	return matr;
}

double* CreateVector(int size) {
	int i;
	double* vect = (double*)malloc(sizeof(double) * size);

	for (i = 0; i < size; i++) {
		vect[i] = (double)(rand() % 10);
	}

	return vect;
}

double GetCond(double* matr, int size) {
	return NormMatrix(matr, size) * NormMatrix(InvMatrix(matr, size), size);
}

double NormMatrix(double* matr, int size) {
	int i, j;
	double norm = 0, tmp = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++)
			tmp += matr[i * size + j] * matr[i * size + j];
	}
	norm = sqrt(tmp);
	return norm;
}

double* InvMatrix(double* matr, int size) {
	int i, j, k;
	double temp;
	double* B = (double*)malloc(sizeof(double) * size * size), * A = (double*)malloc(sizeof(double) * size * size);

	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			A[i * size + j] = matr[i * size + j];

	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++) {
			B[i * size + j] = 0.0;
			if (i == j)
				B[i * size + j] = 1.0;
		}

	for (k = 0; k < size; k++) {
		temp = A[k * size + k];

		for (j = 0; j < size; j++) {
			A[k * size + j] /= temp;
			B[k * size + j] /= temp;
		}

		for (i = k + 1; i < size; i++)
		{
			temp = A[i * size + k];

			for (j = 0; j < size; j++)
			{
				A[i * size + j] -= A[k * size + j] * temp;
				B[i * size + j] -= B[k * size + j] * temp;
			}
		}
	}

	for (k = size - 1; k > 0; k--)
	{
		for (i = k - 1; i >= 0; i--)
		{
			temp = A[i * size + k];

			for (j = 0; j < size; j++)
			{
				A[i * size + j] -= A[k * size + j] * temp;
				B[i * size + j] -= B[k * size + j] * temp;
			}
		}
	}

	return B;
}

double* MatrixMul(double* A, double* X, int size) {
	double* res = (double*)malloc(sizeof(double) * size);
	int i, j;
	for (i = 0; i < size; i++) {
		res[i] = 0;
		for (j = 0; j < size; j++)
			res[i] += A[i * size + j] * X[j];
	}
	return res;
}

void CreateLUMatrix(double* matr, int size, double** L, double** U) {
	double sum;
	int i, j, k;

	for (j = 0; j < size; j++) {
		(*U)[j] = matr[j];
		(*L)[j * size] = matr[j * size] / (*U)[0];
		luOps2++;
	}

	for (i = 1; i < size; i++) {
		for (j = i; j < size; j++) {
			for (k = 0, sum = 0; k < i; k++) {
				sum += (*L)[i * size + k] * (*U)[k * size + j];
				luOps1++;
				luOps2++;
			}
			(*U)[i * size + j] = matr[i * size + j] - sum;
			luOps1++;
			for (k = 0, sum = 0; k < i; k++) {
				sum += (*L)[j * size + k] * (*U)[k * size + i];
				luOps1++;
				luOps2++;
			}
			(*L)[j * size + i] = (matr[j * size + i] - sum) / (*U)[i * size + i];
			luOps2++;
		}
	}
}

double* Solution(double* L, double* U, int size, double* b) {
	double* y = (double*)malloc(sizeof(double) * size), * x = (double*)malloc(sizeof(double) * size), sum;
	int i, j, n = SIZE - 1;

	for (i = 0; i < size; i++) {
		for (j = 0, sum = 0; j < i; j++) {
			sum += L[i * size + j] * y[j];
			luOps1++;
			luOps2++;
		}
		y[i] = b[i] - sum;
		luOps1++;
	}

	for (i = n; i >= 0; i--) {
		for (j = i + 1, sum = 0; j < size; j++) {
			sum += U[i * size + j] * x[j];
			luOps1++;
			luOps2++;
		}
		x[i] = (y[i] - sum) / U[i * size + i];
		luOps1++;
		luOps2++;
	}

	free(y);
	return x;
}

double* CreateAlphaMatrix(double* A, int size) {
	double* res = (double*)malloc(sizeof(double) * size * size);
	int i, j;
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			if (j == i)
				res[i * size + j] = 0;
			else {
				res[i * size + j] = A[i * size + j] / A[i * size + i] * (-1);
				itOps2++;
			}
	return res;
}

double* CreateBetaVector(double* b, double* A, int size) {
	double* res = (double*)malloc(sizeof(double) * size);
	int i;
	for (i = 0; i < size; i++) {
		res[i] = b[i] / A[i * size + i];
		itOps2++;
	}

	return res;
}

double* CreateXVector(double* beta, double* alpha, double* x, int size) {
	double* res = (double*)malloc(sizeof(double) * size);
	int i, j;
	for (i = 0; i < size; i++) {
		res[i] = beta[i];
		for (j = 0; j < size; j++) {
			res[i] += alpha[i * size + j] * x[j];
			itOps1++;
			itOps2++;
		}
	}
	return res;
}

double FindNorm(double* vect, int size) {
	int i;
	double max = 0;
	for (i = 0; i < size; i++)
		if (vect[i] > max)
			max = vect[i];
	return max;
}

double* VectSub(double* v1, double* v2, int size) {
	int i;
	double* res = (double*)malloc(sizeof(double) * size);
	for (i = 0; i < size; i++)
		res[i] = fabs(v1[i] - v2[i]);
	return res;
}

void PrintMatrix(FILE* F, double* matr, int size) {
	int i, j;

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++)
			fprintf(F, "%lf ", matr[i * size + j]);
		fprintf(F, "\n");
	}
	fprintf(F, "\0");

	fseek(F, 0, SEEK_SET);
}

int main(void) {
	double* A, * b, * X, * x;
	double* L, * U;
	double* alpha, * beta, epsilon = 1e-14;
	unsigned int i, tmp, j;
	double cond;
	FILE* F;

	F = fopen("../output1.txt", "w");
	for (j = 0; j < 3000; j++) {

		A = CreateMatrix(SIZE);
		X = CreateVector(SIZE);
		b = MatrixMul(A, X, SIZE);
		cond = GetCond(A, SIZE);
		if (cond <= 0 || cond >= 10000)
			continue;
		fprintf(F, "%.30lf ", cond);

		//--------------------------------------------------------------------------------

		luOps1 = 0;
		luOps2 = 0;

		L = (double*)calloc(SIZE * SIZE, sizeof(double));
		U = (double*)calloc(SIZE * SIZE, sizeof(double));
		CreateLUMatrix(A, SIZE, &L, &U);
		x = Solution(L, U, SIZE, b);
		free(L);
		free(U);
		fprintf(F, "%.15lf ", FindNorm(VectSub(X, x, SIZE), SIZE));
		free(x);

		//--------------------------------------------------------------------------------

		itOps1 = 0;
		itOps2 = 0;

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = calloc(SIZE, sizeof(double));
		for (i = 0; i < SIZE; i++)
			x[i] = beta[i];

		for (i = 0; itOps2 <= luOps2; i++) {
			x = CreateXVector(beta, alpha, x, SIZE);
		}

		fprintf(F, "%.15lf ", FindNorm(VectSub(X, x, SIZE), SIZE));

		free(alpha);
		free(beta);

		//-------------------------------------------------------------------------------

		itOps1 = 0;
		itOps2 = 0;

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = calloc(SIZE, sizeof(double));
		for (i = 0; i < SIZE; i++)
			x[i] = beta[i];

		for (i = 0; itOps2 <= luOps2 / 2; i++) {
			x = CreateXVector(beta, alpha, x, SIZE);
		}

		fprintf(F, "%.15lf ", FindNorm(VectSub(X, x, SIZE), SIZE));

		free(alpha);
		free(beta);

		//-------------------------------------------------------------------------------


		itOps1 = 0;
		itOps2 = 0;

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = calloc(SIZE, sizeof(double));
		for (i = 0; i < SIZE; i++)
			x[i] = beta[i];

		for (i = 0; itOps2 <= luOps2 * 2; i++) {
			x = CreateXVector(beta, alpha, x, SIZE);
		}

		fprintf(F, "%.15lf\n", FindNorm(VectSub(X, x, SIZE), SIZE));

		free(alpha);
		free(beta);

		//-------------------------------------------------------------------------------

		free(A);
		free(b);
		free(X);
		free(x);

	}
	fclose(F);

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	F = fopen("../output2.txt", "w");
	for (epsilon = 1e-13; epsilon < 1; epsilon *= 10) {

		A = CreateMatrix(SIZE);
		X = CreateVector(SIZE);
		b = MatrixMul(A, X, SIZE);
		cond = GetCond(A, SIZE);
		fprintf(F, "%.13lf ", epsilon);

		//--------------------------------------------------------------------------------

		luOps1 = 0;
		luOps2 = 0;

		L = (double*)calloc(SIZE * SIZE, sizeof(double));
		U = (double*)calloc(SIZE * SIZE, sizeof(double));
		CreateLUMatrix(A, SIZE, &L, &U);
		x = Solution(L, U, SIZE, b);
		free(L);
		free(U);
		free(x);
		fprintf(F, "%i %i ", luOps1, luOps2);

		//--------------------------------------------------------------------------------

		itOps1 = 0;
		itOps2 = 0;

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = calloc(SIZE, sizeof(double));
		for (i = 0; i < SIZE; i++)
			x[i] = beta[i];

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++) {
			x = CreateXVector(beta, alpha, x, SIZE);
		}

		fprintf(F, "%i %i\n", itOps1, itOps2);

		free(alpha);
		free(beta);

		//-------------------------------------------------------------------------------

		free(A);
		free(b);
		free(X);
		free(x);
	}
	fclose(F);

	return 0;
}