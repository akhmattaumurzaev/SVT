#include <iostream>
#include <cmath>

#define PI 3.1415926535897932385
#define MAX_N 1e+5

double u_accurate(double x) {
    return std::sin(PI * std::cos(x));
}

double f(double x) {
    return PI * (PI * std::sin(x) * std::sin(x) * std::sin(PI * std::cos(x)) + std::cos(x) * std::cos(PI * std::cos(x)));
}

double* solution(double f(double), int N, double a, double b) {
    double *u = new double[N + 1];
    u[0] = a;
    u[N] = b;

    double *alpha = new double[N - 1], *beta = new double[N - 1];
    double h = 1 / (double) N;
    
    alpha[0] = 1 / 2;
    beta[0] = (f(h) * h * h + u[0]) / 2;

    for (int i = 1; i < N - 1; i++) {
        double temp = 2 - alpha[i - 1];
        alpha[i] = 1 / temp;
        beta[i] = (f((i + 1) * h) / (N * N) + beta[i - 1]) / temp;
    }

    u[N - 1] = (f((N - 1) * h) / (N * N) + b + beta[N - 2]) / (2 - alpha[N - 2]);
    for (int i = N - 2; i >= 0; i--) {
        u[i + 1] = alpha[i] * u[i + 2] + beta[i];
    }

    return u;
}


double C_norm(double u_accurate(double), double *u_calculated, int N) {
    double h = 1 / (double) N;
    double norm = std::abs(u_accurate(h) - u_calculated[1]);
    for (int i = 2; i < N; i++) {
        double temp = std::abs(u_accurate(i * h) - u_calculated[i]);
        if (temp < norm) {
            norm = temp;
        }
    }
    return norm;
}


double L_2_norm(double u_accurate(double), double* u_calculated, int N) {
    double h = 1 / (double) N;
    double norm = 0;
    for (int i = 1; i < N - 1; i++) {
        norm += std::pow(((u_accurate(i * h) - u_calculated[i]) * h), 2);
    }
    return std::sqrt(norm);
}

int main()
{
    for (int N = 10; N < MAX_N + 1; N *= 4) {
        double* u = solution(f, N, u_accurate(0), u_accurate(1));
        std::cout << "N = " << N << std::endl;
        std::cout << "C_norm = " << C_norm(u_accurate, u, N) << std::endl << "L_2_norm = " << L_2_norm(u_accurate, u, N) << std::endl;
        delete u;
    }
}