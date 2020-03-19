#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#pragma warning (disable : 4996)
#define SIZE 12

void PrintMatrix(FILE* F, double* matr, int size);
double* CreateMatrix(int size);
double* MatrixMul(double* A, double* x, int size); //matrix * vector
double* CreateVector(int size);
double FindNorm(double* vect, int size);
double* VectSub(double* v1, double* v2, int size);
double* ReadMatrix(FILE* F);
void PrintVector(FILE* F, double* vect, int size);
double NormMatrix(double* matr, int size);
double* InvMatrix(double* matr, int size);
double GetCond(double* matr, int size);

double* CreateMatrix(int size) {
	int i, j, a = 0;
	double* matr = (double*)malloc(sizeof(double) * size * size);
	//srand((unsigned int)time(0));

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			while ((a = rand() % 10) == 0);
			matr[i * size + j] = (double)a;
		}
	}

	//printf("Matrix with size of %i successfully created\n", size);
	return matr;
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

void PrintVector(FILE* F, double* vect, int size) {
	int i;

	for (i = 0; i < size; i++) {
		fprintf(F, "%lf\n", vect[i]);
	}
	fprintf(F, "\0");

	fseek(F, 0, SEEK_SET);
	//printf("Matrix with size of %i successfully printed\n", size);
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

double* MatrixMul(double* A, double* x, int size) {
	double* res = (double*)malloc(sizeof(double) * size);
	int i, j;

	for (i = 0; i < size; i++) {
		res[i] = 0;
		for (j = 0; j < size; j++)
			res[i] += A[i * size + j] * x[j];
	}

	return res;
}

void CreateLUMatrix(double* matr, int size, double** L, double** U) {
	double sum;
	int i, j, k;

	for (j = 0; j < size; j++) {
		(*U)[j] = matr[j];
		(*L)[j * size] = matr[j * size] / (*U)[0];
	}

	for (i = 1; i < size; i++) {
		for (j = i; j < size; j++) {
			for (k = 0, sum = 0; k < i; k++)
				sum += (*L)[i * size + k] * (*U)[k * size + j];
			(*U)[i * size + j] = matr[i * size + j] - sum;
			for (k = 0, sum = 0; k < i; k++)
				sum += (*L)[j * size + k] * (*U)[k * size + i];
			(*L)[j * size + i] = (matr[j * size + i] - sum) / (*U)[i * size + i];
		}
	}
}

double* Solution(double* L, double* U, int size, double* b) {
	double* y = (double*)malloc(sizeof(double) * size), * x = (double*)malloc(sizeof(double) * size), sum;
	int i, j, n = SIZE - 1;

	for (i = 0; i < size; i++) {
		for (j = 0, sum = 0; j < i; j++) {
			sum += L[i * size + j] * y[j];
		}
		y[i] = b[i] - sum;
	}

	for (i = n; i >= 0; i--) {
		for (j = i + 1, sum = 0; j < size; j++) {
			sum += U[i * size + j] * x[j];
		}
		x[i] = (y[i] - sum) / U[i * size + i];
	}

	//printf("Solution successfully found\n");
	free(y);
	return x;
}

double* CreateVector(int size) {
	int i;
	double* vect = (double*)malloc(sizeof(double) * size);

	for (i = 0; i < size; i++) {
		vect[i] = (double)(rand() % 10);
	}

	//printf("Vector with size of %i successfully created\n", size);
	return vect;
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

double* ChangeVector(double* vect, int size, double num) {
	double* res = (double*)calloc(size, sizeof(double));
	int i;
	for (i = 0; i < size; i++)
		res[i] = vect[i] * (1 + (double)((double)rand() / (double)(RAND_MAX / num)));
	return res;
}

int main(void) {
	double* A, * L, * U, * b, * X, * x, cond, j, * b2;
	int i;
	FILE* F, *F_OUT;
	F = fopen("../../Graphics.txt", "w");
	//F = stdout;

	for (i = 0; i < 10000; i++) {
		L = (double*)calloc(SIZE * SIZE, sizeof(double));
		U = (double*)calloc(SIZE * SIZE, sizeof(double));
		A = CreateMatrix(SIZE);
		X = CreateVector(SIZE);
		b = MatrixMul(A, X, SIZE);
		CreateLUMatrix(A, SIZE, &L, &U);
		x = Solution(L, U, SIZE, b);
		cond = GetCond(A, SIZE);
		if (cond > 0 && cond < 10000) {
			fprintf(F, "%.30lf ", cond);
			fprintf(F, "%.30lf ", FindNorm(VectSub(X, x, SIZE), SIZE));
			fprintf(F, "%.30lf\n", FindNorm(VectSub(MatrixMul(A, x, SIZE), b, SIZE), SIZE));
		}


		free(A);
		free(L);
		free(U);
		free(X);
		free(x);
		free(b);
	}

	fclose(F);
	
	F = fopen("../../Graphics2.txt", "w");

	i = 0;
	do {
		srand(i++);
		A = CreateMatrix(SIZE);
	} while ((cond = GetCond(A, SIZE)) < 1500 && (cond < 0));
	fprintf(stdout, "%.15lf ", cond);
	X = CreateVector(SIZE);

	for (i = 0; i < 1000; i++) {
		L = (double*)calloc(SIZE * SIZE, sizeof(double));
		U = (double*)calloc(SIZE * SIZE, sizeof(double));
		b = MatrixMul(A, X, SIZE);
		j = (double)(rand() % 5 + 1);
		b2 = ChangeVector(b, SIZE, j);
		CreateLUMatrix(A, SIZE, &L, &U);
		x = Solution(L, U, SIZE, b2);

		fprintf(F, "%.15lf ", FindNorm(VectSub(X, x, SIZE), SIZE) / FindNorm(X, SIZE));
		fprintf(F, "%.15f\n", FindNorm(VectSub(b, b2, SIZE), SIZE) / FindNorm(b, SIZE));

		free(L);
		free(U);
		free(x);
	}

	free(A);
	free(b);
	free(X);
	

	fclose(F);
	return 0;
}