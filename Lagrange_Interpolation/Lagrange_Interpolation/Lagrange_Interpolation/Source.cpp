#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

unsigned N;
double Left;
double Right;

// y(x) = 1/(1+25x^2), xi = -1, -0.9 ... 0.9, 1.0; i = 21
double CountFunc(double x) {
    return 1.0 / (double)(1 + 25 * x * x);
}

double CountFuncAngled(double x) {
    double point = 0.3;
    double y = fabs(1.0 / (double)(1 + 25 * x * x) - 1.0 / (double)(1 + 25 * point * point)) + 1.0 / (double)(1 + 25 * point * point);
    return y;
}

void FillChebGrid(vector < double >* xi, vector <double> *yi) {
    for (unsigned i = 0; i <= N; i++)
        (*xi).push_back((Right + Left) / 2 + 0.5 * (Right - Left) * cos(M_PI * (2 * (double)i + 1) / (2 * ((double)N + 1))));
    sort((*xi).begin(), (*xi).end());
    for (unsigned i = 0; i <= N; i++)
        (*yi).push_back(CountFunc((*xi)[i]));
}

void FillChebGridAngled(vector < double >* xi, vector <double>* yi) {
    for (unsigned i = 0; i <= N; i++)
        (*xi).push_back((Right + Left) / 2 + 0.5 * (Right - Left) * cos(M_PI * (2 * (double)i + 1) / (2 * ((double)N + 1))));
    sort((*xi).begin(), (*xi).end());
    for (unsigned i = 0; i <= N; i++)
        (*yi).push_back(CountFuncAngled((*xi)[i]));
}

void FillRandGrid(vector < double >* xi, vector <double>* yi) {
    double h = (Right - Left) / N;
    for (unsigned i = 0; i <= N; i++)
        if (i <= N / 3 || i >= N * 2 / 3)
            (*xi).push_back(Left + i * h);
        else
            (*xi).push_back(Left + i * pow(h, 3));
    sort((*xi).begin(), (*xi).end());
    for (unsigned i = 0; i <= N; i++)
        (*yi).push_back(CountFunc((*xi)[i]));
}

void FillRandGridAngled(vector < double >* xi, vector <double>* yi) {
    double h = (Right - Left) / N;
    for (unsigned i = 0; i <= N; i++)
        if (i <= N / 3 || i >= N * 2 / 3)
            (*xi).push_back(Left + i * h);
        else
            (*xi).push_back(Left + i * pow(h, 3));
    sort((*xi).begin(), (*xi).end());
    for (unsigned i = 0; i <= N; i++)
        (*yi).push_back(CountFuncAngled((*xi)[i]));
}

void FillXTest(vector <double>* x, vector <double> xi) {
    for (unsigned i = 0; i < N; i++)
        (*x).push_back((xi[i + 1] + xi[i]) / 2.0);
}

void FillFuncValues(vector <double>* y, vector <double> x) {
    for (unsigned i = 0; i < N; i++)
        (*y).push_back(CountFunc(x[i]));
}

void FillFuncValuesAngled(vector <double>* y, vector <double> x) {
    for (unsigned i = 0; i < N; i++)
        (*y).push_back(CountFuncAngled(x[i]));
}

void FillPolynomValues(vector <double>* p, vector <double> x, vector <double> xi, vector <double> yi) {
    for (unsigned i = 0; i < N; i++) {  //x loop
        double Lx = 0;
        for (unsigned j = 0; j <= N; j++) { //coeff loop
            double li = 1;
            for (unsigned k = 0; k <= N; k++) { //Li loop
                if (k != j)
                    li = li * (x[i] - xi[k]) / (xi[j] - xi[k]);
            }
            li = li * yi[j];
            Lx += li;
        }
        (*p).push_back(Lx);
    }
}

//Find the biggest difference in values between polynom and function points
double FindMistake(vector <double> func, vector <double> poly) {
    double max = fabs(func[0] - poly[0]);
    for (unsigned i = 1; i < N; i++)
        if (fabs(func[i] - poly[i]) > max)
            max = fabs(func[i] - poly[i]);
    return max;
}

void Processing(string filename, vector <double> xi, vector <double> yi) {

    ofstream file;
    file.open(filename, ios_base::app);

    //filling array of test x values
    vector <double> x;
    FillXTest(&x, xi);

    //filling array of func values in x points
    vector <double> y;
    FillFuncValues(&y, x);

    //filling array of polinom values in x points
    vector <double> p;
    FillPolynomValues(&p, x, xi, yi);

    file << N << " " << FindMistake(y, p) << endl;

    file.close();
}

void ProcessingAngled(string filename, vector <double> xi, vector <double> yi) {

    ofstream file;
    file.open(filename, ios_base::app);

    //filling array of test x values
    vector <double> x;
    FillXTest(&x, xi);

    //filling array of func values in x points
    vector <double> y;
    FillFuncValuesAngled(&y, x);

    //filling array of polinom values in x points
    vector <double> p;
    FillPolynomValues(&p, x, xi, yi);

    file << N << " " << FindMistake(y, p) << endl;

    file.close();
}

int main(void) {

    Left = -1;  //left border
    Right = 1;  //right border

    unsigned MAX = 50;

    //Chebyshev grid
    string filename = "../../output1.txt";

    ofstream file;
    file.open(filename, ios_base::trunc);
    file.close();

    for (N = 1; N < MAX; N++) {

        //filling array of x(i) and y(i) values for chebyshev grid
        vector < double > xi;
        vector < double > yi;
        FillChebGrid(&xi, &yi);

        Processing(filename, xi, yi);
    }
    //------------------------------------------------------------------------------------------

    //Chebyshev grid for angled function
    filename = "../../output2.txt";

    file.open(filename, ios_base::trunc);
    file.close();

    for (N = 1; N < MAX; N++) {

        //filling array of x(i) and y(i) values for chebyshev grid
        vector < double > xi;
        vector < double > yi;
        FillChebGridAngled(&xi, &yi);

        ProcessingAngled(filename, xi, yi);
    }
    //------------------------------------------------------------------------------------------

    //Random grid
    filename = "../../output3.txt";

    file.open(filename, ios_base::trunc);
    file.close();

    for (N = 1; N < MAX; N++) {

        //filling array of x(i) and y(i) values for random grid
        vector < double > xi;
        vector < double > yi;
        FillRandGrid(&xi, &yi);

        Processing(filename, xi, yi);
    }
    //------------------------------------------------------------------------------------------

    //Random grid for angled function
    filename = "../../output4.txt";

    file.open(filename, ios_base::trunc);
    file.close();

    for (N = 1; N < MAX; N++) {

        //filling array of x(i) and y(i) values for random grid
        vector < double > xi;
        vector < double > yi;
        FillRandGridAngled(&xi, &yi);

        ProcessingAngled(filename, xi, yi);
    }

    return 0;
}