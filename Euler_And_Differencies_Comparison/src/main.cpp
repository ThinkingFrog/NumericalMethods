#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

namespace euler_cauchy {
    //xy''+y'+2y=2lnx
    //y''+(1/x)y'+(2/x)y=(2ln(x)/x)
    double LEFT = 1;
    double RIGHT = 2;
    double YLEFT = 0;
    double DEVARLEFT = 1;
    
    double const YLEFT_CONST = YLEFT;
    double const DEVARLEFT_CONST = DEVARLEFT;

    unsigned N;
    double h;

    unsigned imax;
    unsigned coeffCall;
    unsigned rightPartCall;

    //y = ln(x)
    double CountPreciseSolution(double x) {
        double y = log(x);
        return y;
    }

    //y''=2ln(x)/x-z/x-2y/x
    //z = y'
    double CountSecondDevar(double x, double y, double z) {
        double f;
        f = 2 * log(x) / x - z / x - 2 * y / x;
        return f;
    }

    void FillGrid(vector <double>& xi) {
        for (unsigned i = 0; i < N; i++)
            xi.push_back(LEFT + i * h);
    }

    double CountNextY(double x, double y, double z, double h) {
        double yNext;
        yNext =  y + h / 2.0 * (z + z + h * CountSecondDevar(x, y, z));
        return yNext;
    }

    double CountNextZ(double x, double y, double z, double h) {
        double zNext, f = CountSecondDevar(x, y, z);
        zNext = z + h / 2.0 * (f + CountSecondDevar(x + h, y + h * z, z + h * f));
        return zNext;
    }

    void FindSolution(vector <double> xi, vector <double>& yi) {
        vector <double> zi;
        yi.push_back(YLEFT);
        zi.push_back(DEVARLEFT);

        for (unsigned i = 1; i < N; i++) {
            double z = CountNextZ(xi[i - 1], yi[i - 1], zi[i - 1], h);
            rightPartCall += 3;
            coeffCall += 9;
            double y = CountNextY(xi[i - 1], yi[i - 1], zi[i - 1], h);
            rightPartCall++;
            coeffCall += 3;

            zi.push_back(z);
            yi.push_back(y);
        }
    }

    double FindMaxMistake(vector <double> y, vector <double> yi) {
        double max = 0.0;
        imax = 0;
        for (unsigned i = 0; i < y.size(); ++i) {
            double mist = fabs(y[i] - yi[i]);
            if (max < mist) {
                max = mist;
                imax = i;
            }
        }
        return max;
    }
}

namespace finite_differencies {
    //xy''+y'+2y=2lnx
    //y''+(1/x)y'+(2/x)y=(2ln(x)/x)

    const double LEFT = 1;
    const double RIGHT = 2;
    const double A = 0;
    const double B = 1;
    const double C = 1;
    const double D = 0;
    const double alpha = 0;
    const double beta = 0.5;

    unsigned N;
    double h;

    unsigned imax;
    unsigned coeffCall;
    unsigned rightPartCall;

    double error1 = 0;
    double error2 = 0;

    //y = ln(x)
    double CountPreciseSolution(double x) {
        double y = log(x);
        return y;
    }

    void FillXGrid(vector<double>& xi) {
        for (unsigned i = 0; i < N; ++i)
            xi.push_back(LEFT + i * h);
    }

    void FillSupportingVectors(vector <double> xi, vector <double>& pi, vector <double>& qi, vector <double>& fi) {
        for (auto x : xi) {
            pi.push_back(1 / x);
            qi.push_back(2 / x);
            fi.push_back(2 * log(x) / x);

            rightPartCall += 3;
        }
    }

    void FillMatrixVectors(vector <double> pi, vector <double> qi, vector <double> fi, vector <double>& bi, vector <double>& ci, vector <double>& di, vector <double>& ri) {
        for (unsigned i = 0; i < N; i++)
            if (i == 0) {
                bi.push_back(0);
                ci.push_back(h * B - A);
                di.push_back(A);
                ri.push_back(h * (alpha + error1));

                coeffCall += 3;
            }
            else if (i == N - 1) {
                bi.push_back(-ci[i - 1] - 4 * bi[i - 1]);
                ci.push_back(3 * bi[i - 1] - di[i - 1]);
                di.push_back(0);
                ri.push_back(2 * h * bi[i - 1] * (beta + error2) - ri[i - 1]);

                coeffCall += 3;
            }
            else {
                bi.push_back(1 - h / 2 * pi[i]);
                ci.push_back(h * h * qi[i] - 2);
                di.push_back(1 + h / 2 * pi[i]);
                ri.push_back(h * h * fi[i]);
            }
    }

    void FillDeltaVector(vector <double> bi, vector <double> ci, vector <double> di, vector <double>& delta) {
        for (unsigned i = 0; i < N; ++i)
            if (i == 0)
                delta.push_back(-di[0] / ci[0]);
            else
                delta.push_back(-di[i] / (ci[i] + bi[i] * delta[i - 1]));
    }

    void FillLambdaVector(vector <double> bi, vector <double> ci, vector <double> ri, vector <double> delta, vector <double>& lambda) {
        for (unsigned i = 0; i < N; ++i)
            if (i == 0)
                lambda.push_back(ri[0] / ci[0]);
            else
                lambda.push_back((ri[i] - bi[i] * lambda[i - 1]) / (ci[i] + bi[i] * delta[i - 1]));
    }

    void FindSolution(vector <double> delta, vector <double> lambda, vector <double>& y) {
        for (int i = N - 1; i >= 0; --i)
            if (i == N - 1)
                y.push_back(lambda[N - 1]);
            else
                y.push_back(delta[i] * y[N - i - 2] + lambda[i]);

        reverse(y.begin(), y.end());
    }

    double FindMaxMistake(vector <double> y, vector <double> yi) {
        double max = 0.0;
        imax = 0;
        for (unsigned i = 0; i < y.size(); ++i) {
            double mist = fabs(y[i] - yi[i]);
            if (max < mist) {
                max = mist;
                imax = i;
            }
        }
        return max;
    }
}

void Euler_Cauchy_Processing() {
    using namespace euler_cauchy;

    ofstream file1, file2, file3, file4, file5;
    file1.open("../output1_1.txt", ios::out | ios::trunc);
    file2.open("../output2_1.txt", ios::out | ios::trunc);
    file3.open("../output3_1.txt", ios::out | ios::trunc);
    file4.open("../output4_1.txt", ios::out | ios::trunc);
    file5.open("../output5_1.txt", ios::out | ios::trunc);

    for (N = 10; N < 1e4; N += (int)pow(10, (int)log10(N))) {
        imax = 0;

        h = (RIGHT - LEFT) / (N - 1);

        vector <double> xi;
        FillGrid(xi);

        vector <double> y;
        for (auto x : xi)
            y.push_back(CountPreciseSolution(x));

        vector <double> yi;
        FindSolution(xi, yi);

        file1 << h << " " << FindMaxMistake(y, yi) << endl;
        file2 << h << " " << xi[imax] << endl;
        file3 << h << " " << coeffCall + rightPartCall << endl;
    }

    N = 200;
    h = (RIGHT - LEFT) / N;

    vector <double> xi;
    FillGrid(xi);

    vector <double> y;
    for (unsigned i = 0; i < xi.size(); ++i)
        y.push_back(CountPreciseSolution(xi[i]));

    for (double delta = 0.1; delta > 1e-15; delta /= 10) {
        YLEFT = YLEFT_CONST + delta;
        vector <double> yi;
        FindSolution(xi, yi);

        file4 << delta << " " << FindMaxMistake(y, yi) << endl;
    }

    for (double delta = 0.1; delta > 1e-15; delta /= 10) {
        DEVARLEFT = DEVARLEFT_CONST + delta;
        vector <double> yi;
        FindSolution(xi, yi);

        file5 << delta << " " << FindMaxMistake(y, yi) << endl;
    }

    file1.close();
    file2.close();
    file3.close();
    file4.close();
    file5.close();
}

void Finite_Differencies_Processing() {
    using namespace finite_differencies;

    ofstream file1, file2, file3, file4, file5;
    file1.open("../output1_2.txt", ios::out | ios::trunc);
    file2.open("../output2_2.txt", ios::out | ios::trunc);
    file3.open("../output3_2.txt", ios::out | ios::trunc);
    file4.open("../output4_2.txt", ios::out | ios::trunc);
    file5.open("../output5_2.txt", ios::out | ios::trunc);

    for (N = 10; N < 1e4; N += (int)pow(10, (int)log10(N))) {
        h = (RIGHT - LEFT) / (N - 1);

        vector<double> xi;
        FillXGrid(xi);

        vector <double> y;
        for (auto x : xi)
            y.push_back(CountPreciseSolution(x));

        vector <double> pi, qi, fi;
        FillSupportingVectors(xi, pi, qi, fi);

        vector <double> bi, ci, di, ri;
        FillMatrixVectors(pi, qi, fi, bi, ci, di, ri);

        vector <double> delta;
        FillDeltaVector(bi, ci, di, delta);

        vector <double> lambda;
        FillLambdaVector(bi, ci, ri, delta, lambda);

        vector <double> yi;
        FindSolution(delta, lambda, yi);

        file1 << h << " " << FindMaxMistake(y, yi) << endl;
        file2 << h << " " << xi[imax] << endl;
        file3 << h << " " << coeffCall + rightPartCall << endl;
    }

    N = 200;
    h = (RIGHT - LEFT) / (N - 1);
    error1 = 0.1;
    while (error1 >= 1e-15)
    {
        vector<double> xi;
        FillXGrid(xi);

        vector <double> y;
        for (auto x : xi)
            y.push_back(CountPreciseSolution(x));

        vector <double> pi, qi, fi;
        FillSupportingVectors(xi, pi, qi, fi);

        vector <double> bi, ci, di, ri;
        FillMatrixVectors(pi, qi, fi, bi, ci, di, ri);

        vector <double> delta;
        FillDeltaVector(bi, ci, di, delta);

        vector <double> lambda;
        FillLambdaVector(bi, ci, ri, delta, lambda);

        vector <double> yi;
        FindSolution(delta, lambda, yi);

        file4 << error1 << " " << FindMaxMistake(y, yi) << endl;
        error1 /= 10;
    }
    error1 = 0;
    error2 = 0.1;
    while (error2 >= 1e-15)
    {
        vector<double> xi;
        FillXGrid(xi);

        vector <double> y;
        for (auto x : xi)
            y.push_back(CountPreciseSolution(x));

        vector <double> pi, qi, fi;
        FillSupportingVectors(xi, pi, qi, fi);

        vector <double> bi, ci, di, ri;
        FillMatrixVectors(pi, qi, fi, bi, ci, di, ri);

        vector <double> delta;
        FillDeltaVector(bi, ci, di, delta);

        vector <double> lambda;
        FillLambdaVector(bi, ci, ri, delta, lambda);

        vector <double> yi;
        FindSolution(delta, lambda, yi);

        file5 << error2 << " " << FindMaxMistake(y, yi) << endl;
        error2 /= 10;
    }

    file1.close();
    file2.close();
    file3.close();
    file4.close();
    file5.close();
}

int main() {
    Euler_Cauchy_Processing();
    Finite_Differencies_Processing();

    return 0;
}