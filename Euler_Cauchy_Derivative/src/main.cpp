#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

//(x^2+1)y'-2xy=(x^2+1)^2
//y' = f(x,y) = x^2+1+(2xy)/(x^2+1)
//y = x(x^2+1)
double LEFT = 0;
double RIGHT = 2;
double YLEFT = 0;
double N = 20;
double h = (RIGHT - LEFT) / N;
double eps = 1e-1;

//y = x(x^2+1)
double CountPreciseSolution(double x) {
    double y = x * (x * x + 1);
    return y;
}

double CountDevar(double x, double y) {
    double f;
    f = x * x + 1 + (2* x * y) / (x * x + 1);
    return f;
}

void FillGrid(vector <double>& xi) {
    for (unsigned i = 0; i <= N; i++)
        xi.push_back(LEFT + i * h);
}

double CountNextY(double x, double y, double h) {
    double yNext, f = CountDevar(x, y);
    yNext = y + h / 2.0 * (f + CountDevar(x + h, y + h * f));
    return yNext;
}

unsigned FindSolution(vector <double> xi, vector <double>& yi) {
    yi.push_back(YLEFT);
    double y1, y2;
    unsigned jMax = 0;
    for (unsigned i = 1; i <= N; i++) {
        y1 = 100.0;
        y2 = CountNextY(xi[i - 1], yi[i - 1], h);
        unsigned j;
        for (j = 1; fabs(y1 - y2) / 3 > eps;) {
            y1 = y2;
            j *= 2;
            y2 = yi[i - 1];
            for (unsigned k = 0; k < j; ++k)
                y2 = CountNextY(xi[i - 1] + h * k / j, y2, h / j);
        }
        if (jMax < j)
            jMax = j;
        yi.push_back(y2);
    }
    return jMax;
}

double FindMaxMistake(vector <double> y, vector <double> yi) {
    double max = 0.0;
    for (unsigned i = 0; i < y.size(); ++i) {
        double mist = fabs(y[i] - yi[i]);
        if (max < mist)
            max = mist;
    }
    return max;
}

int main(void) {
    vector <double> xi;
    FillGrid(xi);

    vector <double> y;
    for (auto x : xi)
        y.push_back(CountPreciseSolution(x));

    ofstream file1, file2, file3;
    file1.open("../output1.txt", ios::out | ios::trunc);
    file2.open("../output2.txt", ios::out | ios::trunc);
    file3.open("../output3.txt", ios::out | ios::trunc);

    for (eps = 1e-1; eps > 1e-13; eps /= 10) {
        vector <double> yi;
        unsigned iters = FindSolution(xi, yi);

        file1 << eps << " " << FindMaxMistake(y, yi) << endl;
        file2 << eps << " " << iters << endl;
    }

    eps = 1e-5;
    for (YLEFT = 0.1; YLEFT > 1e-10; YLEFT /= 10) {
        vector <double> yi;
        FindSolution(xi, yi);

        file3 << YLEFT << " " << FindMaxMistake(y, yi) << endl;
    }

    file1.close();
    file2.close();
    file3.close();

    return 0;
}