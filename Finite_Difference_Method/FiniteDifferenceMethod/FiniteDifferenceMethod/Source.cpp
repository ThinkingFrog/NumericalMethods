#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

//xy''+y'+2y=2lnx
//y''+(1/x)y'+(2/x)y=(2ln(x)/x)

const double LEFT = 1;
const double RIGHT = 2;
const double A = 0;
const double B = 1;
const double C = 1;
const double D = 0;
const double alpha = 0;
const double beta = log(2);

unsigned N;
double h;

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
    }
}

void FillMatrixVectors(vector <double> pi, vector <double> qi, vector <double> fi, vector <double>& bi, vector <double>& ci, vector <double>& di, vector <double> &ri) {
    for (unsigned i = 0; i < N; i++)
        if (i == 0) {
            bi.push_back(0);
            ci.push_back(h * B - A);
            di.push_back(A);
            ri.push_back(h * alpha + error1);
        }
        else if (i == N - 1) {
            bi.push_back(-ci[i - 1] - 4 * bi[i - 1]);
            ci.push_back(3 * bi[i - 1] - di[i - 1]);
            di.push_back(0);
            ri.push_back(h * bi[i - 1] - ri[i - 1] + error2);
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
    for (unsigned i = 0; i < y.size(); ++i) {
        double mist = fabs(y[i] - yi[i]);
        if (max < mist)
            max = mist;
    }
    return max;
}

int main() {
    ofstream file1, file2, file3;
    file1.open("../../output1.txt", ios::out | ios::trunc);
    file2.open("../../output2.txt", ios::out | ios::trunc);
    file3.open("../../output3.txt", ios::out | ios::trunc);

    for (N = 2; N <= 1e3; N++) {
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
    }

    N = 1000;
    h = (RIGHT - LEFT) / (N - 1);

    for (error1 = 0.1; error1 > 1e-16; error1 /= 10) {

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

        file2 << error1 << " " << FindMaxMistake(y, yi) << endl;
    }

    error1 = 0;
    
    for (error2 = 0.1; error2 > 1e-16; error2 /= 10) {

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

        file3 << error2 << " " << FindMaxMistake(y, yi) << endl;
    }

    file1.close();
    file2.close();
    file3.close();

    return 0;
}