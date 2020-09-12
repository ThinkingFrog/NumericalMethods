#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#pragma warning (disable : 4996)
#define SIZE 12

int cmpfunc(const void* a, const void* b) {
	return (*(int*)a - *(int*)b);
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

double* ReadVector(FILE* F) {
	int i;
	double* vect = (double*)malloc(sizeof(double) * SIZE);
	for (i = 0; i < SIZE; i++)
			fscanf(F, "%lf", &vect[i]);

	//printf("Vector with size of %i successfully read\n", SIZE);
	return vect;
}

double* CreateMatrix(int size) {
	int i, j, a = 0;
	double* matr = (double*)malloc(sizeof(double) * size * size);
	//srand((unsigned int)time(0));

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			a = rand() % 10;
			matr[i * size + j] = (double)a;
		}
	}

	//printf("Matrix with size of %i successfully created\n", size);
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

double* CreateVectorFixed(int size) {
	int i;
	double* vect = (double*)malloc(sizeof(double) * size);

	for (i = 0; i < size; i++) {
		vect[i] = 1;
	}

	//printf("Vector with size of %i successfully created\n", size);
	return vect;
}

double* MatrixMul(double* A, double* X, int size) {	//matrxi-vector multiplication
	double* res = (double*)malloc(sizeof(double) * size);
	int i, j;
	for (i = 0; i < size; i++) {
		res[i] = 0;
		for (j = 0; j < size; j++)
			res[i] += A[i * size + j] * X[j];
	}
	return res;
}

double* VectSub(double* v1, double* v2, int size) {
	int i;
	double* res = (double*)malloc(sizeof(double) * size);
	for (i = 0; i < size; i++)
		res[i] = fabs(v1[i] - v2[i]);
	return res;
}

double* VectMul(double* v1, double num, int size) {	//vect-number multiplication
	int i;
	double* res = (double*)malloc(sizeof(double) * size);
	for (i = 0; i < size; i++)
		res[i] = v1[i] * num;
	return res;
}

double FindNorm(double* vect, int size) {
	int i;
	double max = 0;
	for (i = 0; i < size; i++)
		if (fabs(vect[i]) > max)
			max = fabs(vect[i]);
	return max;
}

double CountLambda (double* x2, double* x1, int size) {
	double tmp1 = 0, tmp2 = 0;
	int i;
	for (i = 0; i < size; i++)
		tmp1 += x1[i];
	tmp1 /= size;
	for (i = 0; i < size; i++)
		tmp2 += x2[i];
	tmp2 /= size;
	return tmp2 / tmp1;
}

int main(void) {
	FILE* F, *F_M;
	double* A, * x1, * x2, lambda, epsilon, lambdapr, lambdafx, *xFixed, norm;
	int i, j, k;
    char* filename[] = { "../Matrix 0 .txt", "../Matrix 1 .txt", "../Matrix 2 .txt", "../Matrix 3 .txt", "../Matrix 4 .txt", "../Matrix 5 .txt", "../Matrix 6 .txt", "../Matrix 7 .txt",
        "../Matrix 8 .txt", "../Matrix 9 .txt", "../Matrix 10 .txt", "../Matrix 11 .txt", "../Matrix 12 .txt", "../Matrix 13 .txt", "../Matrix 14 .txt", "../Matrix 15 .txt", "../Matrix 16 .txt",
        "../Matrix 17 .txt", "../Matrix 18 .txt", "../Matrix 19 .txt", "../Matrix 20 .txt", "../Matrix 21 .txt", "../Matrix 22 .txt", "../Matrix 23 .txt", "../Matrix 24 .txt", "../Matrix 25 .txt",
    "../Matrix 26 .txt", "../Matrix 27 .txt", "../Matrix 28 .txt", "../Matrix 29 .txt", "../Matrix 30 .txt", "../Matrix 31 .txt", "../Matrix 32 .txt", "../Matrix 33 .txt", "../Matrix 34 .txt",
    "../Matrix 35 .txt", "../Matrix 36 .txt", "../Matrix 37 .txt", "../Matrix 38 .txt", "../Matrix 39 .txt", "../Matrix 40 .txt", "../Matrix 41 .txt", "../Matrix 42 .txt", "../Matrix 43 .txt", 
	"../Matrix 44 .txt", "../Matrix 45 .txt", "../Matrix 46 .txt", "../Matrix 47 .txt", "../Matrix 48 .txt", "../Matrix 49 .txt", "../Matrix 50 .txt", "../Matrix 51 .txt", "../Matrix 52 .txt", 
	"../Matrix 53 .txt", "../Matrix 54 .txt", "../Matrix 55 .txt", "../Matrix 56 .txt", "../Matrix 57 .txt", "../Matrix 58 .txt", "../Matrix 59 .txt", "../Matrix 60 .txt", "../Matrix 61 .txt",
		"../Matrix 62 .txt", "../Matrix 63 .txt", "../Matrix 64 .txt", "../Matrix 65 .txt", "../Matrix 66 .txt", "../Matrix 67 .txt", "../Matrix 68 .txt", "../Matrix 69 .txt", "../Matrix 70 .txt", 
	"../Matrix 71 .txt", "../Matrix 72 .txt", "../Matrix 73 .txt", "../Matrix 74 .txt", "../Matrix 75 .txt", "../Matrix 76 .txt", "../Matrix 77 .txt", "../Matrix 78 .txt", "../Matrix 79 .txt", "../Matrix 80 .txt", };

	F = fopen("../Matrix.txt", "r");
	if (F == NULL)
		return 0;
	A = ReadMatrix(F);
	fclose(F);

	F = fopen("../Vector.txt", "r");
	if (F == NULL)
		return 0;
	xFixed = ReadVector(F);
	fclose(F);
	//qsort(xFixed, SIZE, sizeof(double), cmpfunc);

	F = fopen("../Lambda.txt", "r");
	if (F == NULL)
		return 0;
	fscanf(F, "%lf", &lambdafx);
	fclose(F);

	//-------------------------------------------------------------------------------------------------------------------------------------------
	
	F = fopen("../output1_1.txt", "w");
	epsilon = 1e-15;
	while (epsilon < 1) {
		lambda = 0;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		fprintf(F, "%.15lf %.15lf\n", fabs(lambda - lambdafx), epsilon);
		epsilon *= 10;
	}
	fclose(F);

	//------------------------------------------------------------------------------------------

	F = fopen("../output1_2.txt", "w");
	epsilon = 1e-15;
	while (epsilon < 1) {
		lambda = 0;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		//qsort(x1, SIZE, sizeof(double), cmpfunc);
		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(MatrixMul(A, x1, SIZE), VectMul(x1, lambda, SIZE), SIZE), SIZE), epsilon);
		epsilon *= 10;
	}
	lambda = 0;
	epsilon = 1e-15;
	lambdapr = 1;
	fclose(F);

	//------------------------------------------------------------------------------------------

	F = fopen("../output1_3.txt", "w");
	epsilon = 1e-15;
	while (epsilon < 1) {
		lambda = 0;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		qsort(x1, SIZE, sizeof(double), cmpfunc);
		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(x1, xFixed, SIZE), SIZE), epsilon);
		epsilon *= 10;
	}
	lambda = 0;
	epsilon = 1e-15;
	lambdapr = 1;
	fclose(F);
	free(A);

	//-------------------------------------------------------------------------------------------------------------------------------------------

	F = fopen("../output2_1.txt", "w");

    for (k = 0; k < 81; k++) {

        F_M = fopen(filename[k], "r");
        if (F_M == NULL)
            return 0;
        A = ReadMatrix(F_M);
        fclose(F_M);

        lambda = 0;
        epsilon = 1e-13;
        lambdapr = 1;
        x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
        for (i = 0; i < SIZE; i++)
            x1[i] /= norm;
        for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
            lambdapr = lambda;
            x2 = MatrixMul(A, x1, SIZE);
            lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
            free(x1);
            x1 = x2;
        }
        fprintf(F, "%.15lf %i\n", fabs(lambda - lambdafx), k);

		free(A);
	}

	fclose(F);

	//-------------------------------------------------------------------------------------------------------------------------------------------

	F = fopen("../output2_2.txt", "w");

	for (k = 0; k < 81; k++) {

		F_M = fopen(filename[k], "r");
		if (F_M == NULL)
			return 0;
		A = ReadMatrix(F_M);
		fclose(F_M);

		lambda = 0;
		epsilon = 1e-13;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		qsort(x1, SIZE, sizeof(double), cmpfunc);
		fprintf(F, "%.15lf %i\n", FindNorm(VectSub(MatrixMul(A, x1, SIZE), VectMul(x1, lambda, SIZE), SIZE), SIZE), k);

		free(A);
	}

	fclose(F);

	//-------------------------------------------------------------------------------------------------------------------------------------------

	F = fopen("../output2_3.txt", "w");

	for (k = 0; k < 81; k++) {

		F_M = fopen(filename[k], "r");
		if (F_M == NULL)
			return 0;
		A = ReadMatrix(F_M);
		fclose(F_M);

		lambda = 0;
		epsilon = 1e-13;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		qsort(x1, SIZE, sizeof(double), cmpfunc);
		fprintf(F, "%.15lf %i\n", FindNorm(VectSub(x1, xFixed, SIZE), SIZE), k);

		free(A);
	}

	fclose(F);

	//-------------------------------------------------------------------------------------------------------------------------------------------

	F = fopen("../output3_1.txt", "w");

	for (k = 0; k < 16; k += 1) {

		F_M = fopen("../Matrix.txt", "r");
		if (F_M == NULL)
			return 0;
		A = ReadMatrix(F_M);
		fclose(F_M);

		for (i = 0; i < SIZE; i++)
			A[i * SIZE + i] -= (double)k / 100;

		lambda = 0;
		epsilon = 1e-14;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		fprintf(F, "%.15lf %.15lf\n", fabs(lambda - lambdafx), (double)k / 100);

		free(A);
	}

	fclose(F);

	//-------------------------------------------------------------------------------------------------------------------------------------------

	F = fopen("../output3_2.txt", "w");

	for (k = 0; k < 16; k += 1) {

		F_M = fopen("../Matrix.txt", "r");
		if (F_M == NULL)
			return 0;
		A = ReadMatrix(F_M);
		fclose(F_M);

		for (i = 0; i < SIZE; i++)
			A[i * SIZE + i] -= (double)k / 100;

		lambda = 0;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		qsort(x1, SIZE, sizeof(double), cmpfunc);
		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(MatrixMul(A, x1, SIZE), VectMul(x1, lambda, SIZE), SIZE), SIZE), (double)k / 100);

		free(A);
	}

	fclose(F);

	//-------------------------------------------------------------------------------------------------------------------------------------------

	F = fopen("../output3_3.txt", "w");

	for (k = 0; k < 16; k += 1) {

		F_M = fopen("../Matrix.txt", "r");
		if (F_M == NULL)
			return 0;
		A = ReadMatrix(F_M);
		fclose(F_M);

		for (i = 0; i < SIZE; i++)
			A[i * SIZE + i] -= (double)k / 100;

		lambda = 0;
		lambdapr = 1;
		x1 = CreateVectorFixed(SIZE);
		norm = fabs(FindNorm(x1, SIZE));
		for (i = 0; i < SIZE; i++)
			x1[i] /= norm;
		for (i = 0; fabs(lambda - lambdapr) > epsilon; i++) {
			lambdapr = lambda;
			x2 = MatrixMul(A, x1, SIZE);
			lambda = CountLambda(x2, x1, SIZE);
			norm = fabs(FindNorm(x2, SIZE));
			for (j = 0; j < SIZE; j++)
				x2[j] /= norm;
			free(x1);
			x1 = x2;
		}
		qsort(x1, SIZE, sizeof(double), cmpfunc);
		fprintf(F, "%.15lf %.15lf\n", FindNorm(VectSub(x1, xFixed, SIZE), SIZE), (double)k / 100);

		free(A);
	}

	fclose(F);

	return 0;
}