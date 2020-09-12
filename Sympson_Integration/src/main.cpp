#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>

using namespace std;

unsigned N = 2;
double Left = -1;
double Right = 1;
double RightFixed = Right;
double LeftFixed = Left;
double eps;
double integralFixed = -5.333333333333333;
double DevarMin = -262.8570;
double DevarMax = 2760;

double CountFunc(double x);
void FillXiVector(vector < double >* xi, double h);
double CountIntegral(void);
double CountFuncDevar(double x);

//filling array of x(i) values
void FillXiVector(vector < double >* xi, double h) {
    (*xi).clear();
    for (unsigned i = 0; i < N + 1; i++) {
        double x = Left + i * h;
        (*xi).push_back(x);
    }
}

// y(x) = 7x^6 - 4x^5 - 10x^4 + 17x^3 + x^2 + 21x - 2
double CountFunc(double x) {
    double y = 7 * pow(x, 6) - 4 * pow(x, 5) - 10 * pow(x, 4) + 17 * pow(x, 3) + x * x + 21 * x - 2;
    return y;
}

// y(IV)(x) = 2520 * x^2 - 480 * x - 240
double CountFuncDevar(double x) {
    double y = 2520 * x * x - 480 * x - 240;
    return y;
}

double CountIntegral(void) {
    double res = 0;
    double h = (Right - Left) / N;
    double B = 1.0 / 3.0;
    array < double, 3 > a = { 1, 4, 1 };
    vector < double > xi;

    FillXiVector(&xi, h);

    for (unsigned i = 0; i < N + 1; i++) {
        res += a[i] * CountFunc(xi[i]);
    }
    res *= B * h;

    return res;
}

int main(void) {

    ofstream file1, file2;
    file1.open("../output1_2.txt");
    file2.open("../output2_2.txt");
    file1.precision(15);
    file2.precision(15);

    for (eps = 1e-1; eps > 1e-15; eps /= 10) {

        double integral = 0, integralPrev = 100, step;
        unsigned M, L;

        for (M = 1, L = 0; ; M *= 2, L++) {
            integralPrev = integral;
            integral = 0;

            Right = RightFixed;
            Left = LeftFixed;
            step = (Right - Left) / (double)M;
            Right -= step * (M - 1);

            for (unsigned k = 1; k <= M; k++, Right += step, Left += step) {
                integral += CountIntegral();
            }

            if ((fabs(integral - integralPrev) / (15.0)) < eps)
                break;
        }

        double h = (Right - Left) / N;
        double epsTheorMax = pow(step, 4) / 90.0 * DevarMax;
        double epsTheorMin = -pow(step, 5) / 90.0 * DevarMin;

        file1 << fabs(integralFixed - integral) << " " << eps << " " << epsTheorMax << " " << epsTheorMin << endl;
        file2 << L << " " << eps << endl;

    }

    file1.close();
    file2.close();

    return 0;
}