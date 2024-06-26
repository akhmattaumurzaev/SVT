#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

#define PI 3.1415926535897932385
#define MAX_N 1e+3


double u_accurate(double x, double y) {
    return sin(5.0 * x) * sin(5.0 * y);
}


double f(double x, double y) {
    return -50.0 * std::sin(5.0 * x) * std::sin(5.0 * y);
}


double C_norm(double u_accurate(double, double), Sparse::Vector &u_calculated, int N) {
    double norm = 0;
    for (int i = 0; i < (N - 1) * (N - 1); i++) {
        double temp = std::abs(u_accurate((double)(i % (N - 1) + 1) / N, (double)(i / (N - 1) + 1) / N) - u_calculated[i]);
        if (temp > norm) {
            norm = temp;
        }
    }
    return norm;
}


double L_2_norm(double u_accurate(double, double), Sparse::Vector &u_calculated, int N) {
    double norm = 0;
    for (int i = 0; i < (N - 1) * (N - 1); i++) {
        norm += std::pow((u_accurate((double)(i % (N - 1) + 1) / N, (double)(i / (N - 1) + 1) / N) - u_calculated[i]), 2);
    }
    return std::sqrt(norm) / N;
}

int main(int argc, char *argv[])
{
    for (int N = 10; N < 500; N *= 2) {
        Sparse::Matrix A;
        Sparse::Vector b;
        Sparse::Vector sol;

        A.SetInterval(0, (N - 1) * (N - 1));
        b.SetInterval(0, (N - 1) * (N - 1));
        sol.SetInterval(0, (N - 1) * (N - 1));

        Solver S(Solver::INNER_DDPQILUC);
        S.SetParameter("absolute_tolerance", "1e-10");
        S.SetParameter("relative_tolerance", "1e-6");

        for(int i = 0; i < (N - 1) * (N - 1); i++) {
            A[i][i] = 4;
            b[i] = -f((double) (i % (N - 1) + 1) / N, (double) (i / (N - 1) + 1) / N) / N / N;           

            if (i % (N - 1) + 1 < (N - 1)) {
                A[i][i + 1] = -1;
            } else {
                b[i] += u_accurate(1, (double) (i / (N - 1) + 1) / N);
            }

            if (i % (N - 1) > 0) {
                A[i][i - 1] = -1;
            }

            if (i + N - 1 < (N - 1) * (N - 1)) {
                A[i][i + N - 1] = -1;
            }
            else {
                b[i] += u_accurate((double) (i % (N - 1) + 1) / N, 1);
            }

            if (i - N > 0) {
                A[i][i - N + 1] = -1;
            }
        }

        S.SetMatrix(A);

        bool solved = S.Solve(b, sol);

        std::cout << "N = " << N << std::endl;
        std::cout << "C_norm = " << C_norm(u_accurate, sol, N) << std::endl << "L_2_norm = " << L_2_norm(u_accurate, sol, N) << std::endl;
        cout << "num.iters: " << S.Iterations() << endl;
        cout << "prec.time: " << S.PreconditionerTime() << endl;
        cout << "iter.time: " << S.IterationsTime() << endl;
        if(!solved){
            cout << "Linear solver failure!" << endl;
            cout << "Reason: " << S.ReturnReason() << endl;
        }
    }
	return 0;
}