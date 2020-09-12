#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

unsigned N = 3;
double Left = 1;
double Right = 4;
const double RightFixed = Right;
const double LeftFixed = Left;
double eps;
double integralFixed = 12863.25;

// y(x) = 7x^6 - 4x^5 - 10x^4 + 17x^3 + x^2 + 21x - 2
double CountFunc(double x) {
    double y = 7 * pow(x, 6) - 4 * pow(x, 5) - 10 * pow(x, 4) + 17 * pow(x, 3) + x * x + 21 * x - 2;
    return y;
}

double CountIntegral() {
    double res = 0;
    std::vector<double> ti = {-1.0, -0.289897948556635619639456814941, 0.689897948556635619639456814941};
    std::vector<double> A = {0.222222222222222222222222222222, 1.02497165237684322767762689304, 0.752806125400934550100150884739};

    for (unsigned i = 0; i < N; ++i)
        res += (Right - Left) / 2.0 * A[i] * CountFunc((Left + Right) / 2.0 + (Right - Left) / 2.0 * ti[i]);

    return res;
}

int main() {
    std::ofstream file1, file2, file3;
    
    file1.open("../output1.txt");
    file2.open("../output2.txt");
    file3.open("../output3.txt");

    file1.precision(30);
    file2.precision(30);
    file3.precision(30);

    for (eps = 1e-1; eps > 1e-11; eps /= 10) {
        double integral = 0, integralPrev = 100, step;
        unsigned M, L;

        for (M = 1, L = 0; ; M *= 2, ++L) {
            integralPrev = integral;
            integral = 0;

            Right = RightFixed;
            Left = LeftFixed;
            step = (Right - Left) / (double)M;
            Right = Left + step;

            for (unsigned k = 1; k <= M; k++, Right += step, Left += step)
                integral += CountIntegral();

            if ((fabs(integral - integralPrev) / (31.0)) < eps)
                break;
        }

        file1 << eps << "\t" << fabs(integral - integralFixed) << std::endl;
        file2 << eps << "\t" << L << std::endl;
        file3 << eps << "\t" << fabs(integral - integralFixed) << "\t" << pow(step, 5) << std::endl;
    }

    return 0;
}