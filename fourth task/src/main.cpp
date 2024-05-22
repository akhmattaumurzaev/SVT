#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;

double f(double x, double y) {
    return sin(5.0 * x) * sin(5.0 * y);
}

int main()
{
    for (int N = 128; N < 8193; N *= 2) {
        // Create sparse matrix, RHS vector and solution vector
        Sparse::Matrix A;
        Sparse::Vector b;
        Sparse::Vector sol;
        // Set their size
        A.SetInterval(0, (N - 1) * (N - 1));
        b.SetInterval(0, (N - 1) * (N - 1));
        sol.SetInterval(0, (N - 1) * (N - 1));

        // Get solver
        // All inner INMOST solvers are BiCGStab
        // with different preconditioners, let's use ILU2
        Solver S(Solver::INNER_ILU2);
        S.SetParameter("absolute_tolerance", "1e-10");
        S.SetParameter("relative_tolerance", "1e-6");

        for(int i = 0; i < (N - 1) * (N - 1); i++) {
            A[i][i] = 4;
            b[i] = -50.0 * f((double) (i % (N - 1) + 1) / N, (double) (i / (N - 1) + 1) / N) / N / N;           

            if (i % (N - 1) + 1 < (N - 1)) {
                A[i][i + 1] = -1;
            } else {
                b[i] += f(1, (double) (i / (N - 1) + 1) / N);
            }

            if (i % (N - 1) > 0) {
                A[i][i - 1] = -1;
            }

            if (i + N - 1 < (N - 1) * (N - 1)) {
                A[i][i + N - 1] = -1;
            }
            else {
                b[i] += f((double) (i % (N - 1) + 1) / N, 1);
            }

            if (i - N > 0) {
                A[i][i - N + 1] = -1;
            }
        }

        // Set matrix in the solver;
        // this also computes preconditioner
        S.SetMatrix(A);

        // Solve
        bool solved = S.Solve(b, sol);
        cout << "num.iters: " << S.Iterations() << endl;
        cout << "prec.time: " << S.PreconditionerTime() << endl;
        cout << "iter.time: " << S.IterationsTime() << endl;
        if(!solved){
            cout << "Linear solver failure! N = " << N << endl;
            cout << "Reason: " << S.ReturnReason() << endl;
            return 0;
        }

        double L2_error = 0.0;
        double C_error = 0.0;
        double temp = 0.0;
        for (int i = 0; i < (N - 1) * (N - 1); i++) {
            temp = f((double)(i % (N - 1) + 1) / N, (double)(i / (N - 1) + 1) / N);
            L2_error += pow((sol[i] - temp), 2);
            if (abs(sol[i] - temp) > C_error) {
                C_error = temp;
            }
        }
        cout << "L2 error = " << L2_error << endl;
        cout << "C error = " << C_error << endl;

    }
    return 0;
}