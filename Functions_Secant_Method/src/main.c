#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning (disable : 4996)

typedef struct {
	double x, y;
} point_t;

double SolveFunc(double x, int funcNum) {	//returns solution of designated function with given number x
	switch (funcNum) {
	default:
		printf("Incorrect function number\n");
		return 0;
		break;
	case 1:
		return (3 * x * x * x * x + 4 * x * x * x - 12 * x * x + 1);
		break;
	case 2:
		return (6 * x + 3 - pow(5.0, x));
		break;
	case 3:
		return (2 * x * x - 8);
		break;
	}
}

point_t SearchNextPoint(point_t c1, point_t c2, int funcNum) {	//returns next point for hordes method
	point_t res;
	res.x = (c2.x - (c1.x - c2.x) * SolveFunc(c2.x, funcNum) / (SolveFunc(c1.x, funcNum) - SolveFunc(c2.x, funcNum)));
	res.y = SolveFunc(res.x, funcNum);
	return res;
}

point_t DivisionSearchNextPoint(point_t c1, point_t c2, int funcNum) {	//returns next point for division method
	point_t res;
	res.x = (c2.x + c1.x) / 2;
	res.y = SolveFunc(res.x, funcNum);
	return res;
}

void HordesSolution(point_t p1, point_t p2, double e, int funcNum, FILE* F) {
	int r = 0;
	while (fabs(p2.x - SearchNextPoint(p1, p2, funcNum).x) >= e) {
		p2 = SearchNextPoint(p1, p2, funcNum);
		r++;
	}
	fprintf(F, /*"%.15f, %.15f with r = "*/"%i "/*and e = */"%.15f\n", /*p2.x, p2.y, */r, e);
}

void DivisionSolution(point_t p1, point_t p2, double e, int funcNum, FILE* F) {
	int r = 0;
	point_t p3;
	while (fabs(p2.x - p1.x)/*SearchNextPoint(p1, p2, funcNum).x)*/ >= 2 * e) {
		p3 = DivisionSearchNextPoint(p1, p2, 1);
		if (p3.y * p1.y > 0)
			p1 = p3;
		else
			p2 = p3;
		r++;
	}
	fprintf(F, /*"%.15f, %.15f with r = "*/"%i "/*and e = */"%.15f\n", /*p2.x, p2.y, */r, e);
}

int main(void) {
	FILE* F;
	point_t p1, p2;
	double e;
	int funcNum;
	F = fopen("../output1_1.txt", "w");
	//F = stdout;

	//fprintf(stdout, "First function: 3x^4 + 4x^3 - 12x^2 + 1 = 0\n");
	//fprintf(stdout, "\nHordes method solution:\n");

	//fprintf(stdout, "\nLooking between points (-3, 28) and (-2.5, -19.31) with floating e\n");
	funcNum = 1;
	e = 1e-1;
	p1.x = -3;
	p1.y = SolveFunc(p1.x, 1);
	p2.x = -2.5;
	p2.y = SolveFunc(p2.x, 1);
	while (e >= 1e-15) {
		HordesSolution(p1, p2, e, funcNum, F);
		e /= 10;
	}

	F = fopen("../output1_2.txt", "w");
	//F = stdout;
	//fprintf(F, "\nLooking between floating points with fixed e = 10^-14\n");
	e = 1e-14;
	p1.x = -5;
	p1.y = SolveFunc(p1.x, 1);
	p2.x = 0;
	p2.y = SolveFunc(p2.x, 1);
	while (fabs(p2.x - p1.x) > 0.25) {
		fprintf(F, "%.4f ", fabs(p1.x - p2.x));
		HordesSolution(p1, p2, e, funcNum, F);
		p1.x += 0.125;
		p1.y = SolveFunc(p1.x, 1);
		p2.x -= 0.125;
		p2.y = SolveFunc(p2.x, 1);
	}

	F = fopen("../output1_3.txt", "w");
	//F = stdout;
	//fprintf(F, "\nDivision method solution:\n");

	//fprintf(F, "\nLooking between points (-3, 28) and (-2.5, -19.31) with floating e\n");
	funcNum = 1;
	e = 1e-1;
	p1.x = -3;
	p1.y = SolveFunc(p1.x, 1);
	p2.x = -2.5;
	p2.y = SolveFunc(p2.x, 1);
	while (e >= 1e-15) {
		DivisionSolution(p1, p2, e, funcNum, F);
		e /= 10;
	}

	F = fopen("../output1_4.txt", "w");
	//F = stdout;
	//fprintf(F, "\nLooking between floating points with fixed e = 10^-14\n");
	e = 1e-14;
	p1.x = -5;
	p1.y = SolveFunc(p1.x, 1);
	p2.x = 0;
	p2.y = SolveFunc(p2.x, 1);
	while (fabs(p2.x - p1.x) > 0.25) {
		fprintf(F, "%.2f ", fabs(p1.x - p2.x));
		DivisionSolution(p1, p2, e, funcNum, F);
		p1.x += 0.125;
		p1.y = SolveFunc(p1.x, 1);
		p2.x -= 0.125;
		p2.y = SolveFunc(p2.x, 1);
	}

	F = fopen("../output2_1.txt", "w");
	//F = stdout;
	//fprintf(F, "\n\nSecond function: 5^x = 6x + 3\n");
	//fprintf(F, "\nHordes method solution:\n");

	//fprintf(F, "\nLooking between points (1.5, 0.82) and (2, -10) with floating e\n");
	funcNum = 2;
	e = 1e-1;
	p1.x = 1.5;
	p1.y = SolveFunc(p1.x, 2);
	p2.x = 2;
	p2.y = SolveFunc(p2.x, 2);
	while (e >= 1e-15) {
		HordesSolution(p1, p2, e, funcNum, F);
		e /= 10;
	}

	F = fopen("../output2_2.txt", "w");
	//F = stdout;
	//fprintf(F, "\nLooking between floating points with fixed e = 10^-14\n");
	e = 1e-14;
	p1.x = 1;
	p1.y = SolveFunc(p1.x, 2);
	p2.x = 3;
	p2.y = SolveFunc(p2.x, 2);
	while (fabs(p2.x - p1.x) > 0.25) {
		fprintf(F, "%.2f ", fabs(p1.x - p2.x));
		HordesSolution(p1, p2, e, funcNum, F);
		p1.x += 0.125;
		p1.y = SolveFunc(p1.x, 2);
		p2.x -= 0.125;
		p2.y = SolveFunc(p2.x, 2);
	}

	F = fopen("../output2_3.txt", "w");
	//F = stdout;
	//fprintf(F, "\nDivision method solution:\n");

	//fprintf(F, "\nLooking between points (1.5, 0.82) and (2, -10) with floating e\n");
	funcNum = 2;
	e = 1e-1;
	p1.x = 1.5;
	p1.y = SolveFunc(p1.x, 2);
	p2.x = 2;
	p2.y = SolveFunc(p2.x, 2);
	while (e >= 1e-15) {
		DivisionSolution(p1, p2, e, funcNum, F);
		e /= 10;
	}

	F = fopen("../output2_4.txt", "w");
	//F = stdout;
	//fprintf(F, "\nLooking between floating points with fixed e = 10^-14\n");
	e = 1e-14;
	p1.x = 1;
	p1.y = SolveFunc(p1.x, 2);
	p2.x = 3;
	p2.y = SolveFunc(p2.x, 2);
	while (fabs(p2.x - p1.x) > 0.25) {
		fprintf(F, "%.2f ", fabs(p1.x - p2.x));
		DivisionSolution(p1, p2, e, funcNum, F);
		p1.x += 0.125;
		p1.y = SolveFunc(p1.x, 2);
		p2.x -= 0.125;
		p2.y = SolveFunc(p2.x, 2);
	}

	fclose(F);
	return 0;
}