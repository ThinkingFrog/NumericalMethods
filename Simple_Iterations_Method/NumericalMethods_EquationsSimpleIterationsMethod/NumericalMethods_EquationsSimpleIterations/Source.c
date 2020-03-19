#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#pragma warning (disable : 4996)
#define SIZE 12

void PrintMatrix(FILE* F, double* matr, int size);
void PrintVector(FILE* F, double* vect, int size);
double* CreateMatrix(int size);
double* MatrixMul(double* A, double* X, int size); //matrix * vector
double* CreateVector(int size);
double FindNorm(double* vect, int size);
double* VectSub(double* v1, double* v2, int size);
double* ReadMatrix(FILE* F);
double* CreateBetaVector(double* b, double* A, int size);
double* CreateXVector(double* beta, double* alpha, double* x, int size);
double* CreateAlphaMatrix(double* A, int size);
double NormMatrix(double* matr, int size);
double GetDet(double* A, int size);

double GetDet(double* matr, int size) {
	int i, j, k;
	double sum, mul;
	sum = 0;

	if (size == 2)
		return matr[0] * matr[3] - matr[1] * matr[2];

	for (i = 0; i < size; i++) {
		mul = 1;
		k = i;
		for (j = 0; j < size; j++) {
			mul *= matr[j * size + k];
			k = (k == size - 1 ? 0 : k + 1);
		}
		sum += mul;
	}
	for (i = 0; i < size; i++) {
		mul = 1;
		k = i;
		for (j = 0; j < size; j++) {
			mul *= matr[j * size + k];
			k = (k == 0 ? size - 1 : k - 1);
		}
		sum -= mul;
	}
	return sum;
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

void PrintVector(FILE* F, double* vect, int size) {
	int i;

	for (i = 0; i < size; i++) {
		fprintf(F, "%lf\n", vect[i]);
	}
	fprintf(F, "\0");

	fseek(F, 0, SEEK_SET);
	//printf("Matrix with size of %i successfully printed\n", size);
}

double* CreateMatrix(int size) {
	int i, j, a = 0;
	double* matr = (double*)malloc(sizeof(double) * size * size);
	//srand((unsigned int)time(0));

	for (i = 0; i < size; i++) {
		while ((a = rand() % 10) <= 1);
		matr[i * size + i] = (double)(a);
		for (j = 0; j < size; j++) {
			if (j != i) {
				while ((a = rand() % 10) == 0 || a >= matr[i * size + i]);
				matr[i * size + j] = (double)(a) / 10;
			}
		}
	}

	//printf("Matrix with size of %i successfully created\n", size);
	return matr;
}

double* ReadMatrix(FILE* F) {
	int i, j;
	double* matr = (double*)malloc(sizeof(double) * SIZE * SIZE);
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			fscanf(F, "%lf", &matr[i * SIZE + j]);

	//printf("Matrix with size of %i successfully read\n", SIZE);
	return matr;
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
	//printf("Matrix with size of %i successfully printed\n", size);
}

double* CreateVector(int size) {
	int i;
	double* vect = (double*)malloc(sizeof(double) * size);
	//srand((unsigned int)time(0));

	for (i = 0; i < size; i++) {
		vect[i] = (double)(rand() % 10);
	}

	//printf("Vector with size of %i successfully created\n", size);
	return vect;
}

double* CreateAlphaMatrix(double* A, int size) {
	double* res = (double*)malloc(sizeof(double) * size * size);
	int i, j;
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			if (j == i)
				res[i * size + j] = 0;
			else
				res[i * size + j] = A[i * size + j] / A[i * size + i] * (-1);
	return res;
}

double* CreateBetaVector(double* b, double* A, int size) {
	double* res = (double*)malloc(sizeof(double) * size);
	int i;
	for (i = 0; i < size; i++)
		res[i] = b[i] / A[i * size + i];

	return res;
}

double* CreateXVector(double* beta, double* alpha, double* x, int size) {
	double* res = (double*)malloc(sizeof(double) * size);
	int i, j;
	for (i = 0; i < size; i++) {
		res[i] = beta[i];
		for (j = 0; j < size; j++)
			res[i] += alpha[i * size + j] * x[j];
	}
	return res;
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

int main(void) {
	FILE* F;
	double* A, * b, * alpha, * beta, * x, * X, epsilon = 1e-1, * tmpA, * tmpb, k;
	int i, j;

	A = CreateMatrix(SIZE);
	X = CreateVector(SIZE);
	b = MatrixMul(A, X, SIZE);
	tmpA = (double*)malloc(sizeof(double) * SIZE * SIZE);
	tmpb = (double*)malloc(sizeof(double) * SIZE);
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			tmpA[i * SIZE + j] = A[i * SIZE + j];
	for (i = 0; i < SIZE; i++)
		tmpb[i] = b[i];

	//---------------------------------------------------------------------------------------------

		// ||x* - x|| (eps)
	F = fopen("../../Dependings1.txt", "w");
	//F = stdout;
	while (epsilon >= 1e-15) {

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = beta;

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++)
			x = CreateXVector(beta, alpha, x, SIZE);

		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(X, x, SIZE), SIZE), epsilon);

		free(alpha);
		free(beta);
		free(x);

		epsilon /= 10;
	}
	fclose(F);

	//---------------------------------------------------------------------------------------------


		// ||Ax - b|| (eps)
	epsilon = 1e-1;
	F = fopen("../../Dependings2.txt", "w");
	//F = stdout;
	while (epsilon >= 1e-15) {

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = beta;

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++)
			x = CreateXVector(beta, alpha, x, SIZE);

		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(MatrixMul(A, x, SIZE), b, SIZE), SIZE), epsilon);

		free(alpha);
		free(beta);
		free(x);

		epsilon /= 10;
	}
	fclose(F);

	//---------------------------------------------------------------------------------------------


		// ||x* - x|| (det(A))
	epsilon = 1e-14;
	F = fopen("../../Dependings3.txt", "w");
	//	F = stdout;
	while (GetDet(A, SIZE) > epsilon) {

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = beta;

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++)
			x = CreateXVector(beta, alpha, x, SIZE);

		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(X, x, SIZE), SIZE), GetDet(A, SIZE));

		free(alpha);
		free(beta);
		free(x);

		for (i = 0; i < SIZE; i++)
			for (j = 0; j < SIZE; j++)
				A[i * SIZE + j] /= (double)pow((double)2, (double)1 / (double)SIZE);

		for (i = 0; i < SIZE; i++)
			b[i] /= (double)pow((double)2, (double)1 / (double)SIZE);
	}
	fclose(F);
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			A[i * SIZE + j] = tmpA[i * SIZE + j];
	for (i = 0; i < SIZE; i++)
		b[i] = tmpb[i];

	//---------------------------------------------------------------------------------------------

		// ||Ax - b|| (det(A))
	epsilon = 1e-14;
	F = fopen("../../Dependings4.txt", "w");
	//F = stdout;
	while (GetDet(A, SIZE) > epsilon) {

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = beta;

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++)
			x = CreateXVector(beta, alpha, x, SIZE);

		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(MatrixMul(A, x, SIZE), b, SIZE), SIZE), GetDet(A, SIZE));

		free(alpha);
		free(beta);
		free(x);

		for (i = 0; i < SIZE; i++)
			for (j = 0; j < SIZE; j++)
				A[i * SIZE + j] /= (double)pow((double)2, (double)1 / (double)SIZE);

		for (i = 0; i < SIZE; i++)
			b[i] /= (double)pow((double)2, (double)1 / (double)SIZE);
	}
	fclose(F);
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			A[i * SIZE + j] = tmpA[i * SIZE + j];
	for (i = 0; i < SIZE; i++)
		b[i] = tmpb[i];

	//---------------------------------------------------------------------------------------------

		// N (eps)
	epsilon = 1e-1;
	F = fopen("../../Dependings5.txt", "w");
	//F = stdout;
	while (epsilon >= 1e-15) {

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = beta;

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++)
			x = CreateXVector(beta, alpha, x, SIZE);

		fprintf(F, "%i %.15lf\n", i, epsilon);

		free(alpha);
		free(beta);
		free(x);

		epsilon /= 10;
	}
	fclose(F);

	//---------------------------------------------------------------------------------------------

		// N (x0)
	epsilon = 1e-14;
	k = 100;
	F = fopen("../../Dependings6.txt", "w");
	//F = stdout;
	while (k >= 1) {

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = (double*)malloc(sizeof(double) * SIZE); 
		for (i = 0; i < SIZE; i++)
			x[i] = X[i];
		for (i = 0; i < SIZE; i++)
			x[i] *= k;

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++)
			x = CreateXVector(beta, alpha, x, SIZE);

		fprintf(F, "%i %.15lf\n", i, k);

		free(alpha);
		free(beta);
		free(x);

		k /= 1.1;
	}
	fclose(F);
	//---------------------------------------------------------------------------------------------


		// N (det(A))
	epsilon = 1e-14;
	F = fopen("../../Dependings7.txt", "w");
	//F = stdout;
	while (GetDet(A, SIZE) > epsilon) {

		alpha = CreateAlphaMatrix(A, SIZE);
		beta = CreateBetaVector(b, A, SIZE);
		x = beta;

		for (i = 0; FindNorm(VectSub(X, x, SIZE), SIZE) > epsilon; i++)
			x = CreateXVector(beta, alpha, x, SIZE);

		fprintf(F, "%i %.15lf\n", i, GetDet(A, SIZE));

		free(alpha);
		free(beta);
		free(x);

		for (i = 0; i < SIZE; i++)
			for (j = 0; j < SIZE; j++)
				A[i * SIZE + j] /= (double)pow((double)2, (double)1 / (double)SIZE);

		for (i = 0; i < SIZE; i++)
			b[i] /= (double)pow((double)2, (double)1 / (double)SIZE);
	}
	fclose(F);
	for (i = 0; i < SIZE; i++)
		for (j = 0; j < SIZE; j++)
			A[i * SIZE + j] = tmpA[i * SIZE + j];
	for (i = 0; i < SIZE; i++)
		b[i] = tmpb[i];

	//---------------------------------------------------------------------------------------------
	return 0;
}