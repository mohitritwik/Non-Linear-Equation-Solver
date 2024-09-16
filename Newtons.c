#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double** jacobian_exercise(double x, double y, double z) {
    double** jacobian = (double**)malloc(3 * sizeof(double*));
    for (int i = 0; i < 3; i++) {
        jacobian[i] = (double*)malloc(3 * sizeof(double));
    }
    jacobian[0][0] = 1;
    jacobian[0][1] = 1;
    jacobian[0][2] = 1;
    jacobian[1][0] = 2;
    jacobian[1][1] = 2 * y;
    jacobian[1][2] = 2 * z;
    jacobian[2][0] = exp(x);
    jacobian[2][1] = x;
    jacobian[2][2] = -x;
    return jacobian;
}

double* function_exercise(double x, double y, double z) {
    double* result = (double*)malloc(3 * sizeof(double));
    result[0] = x + y + z - 3;
    result[1] = 2 * x + y * y + 2 * z - 5;
    result[2] = exp(x) + x * y - x * z - 1;
    return result;
}
void forward_elimination(double **A, double *b, int n) {
    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
}

double *backward_substitution(double **A, double *b, int n) {
    double *x = (double *)malloc(n * sizeof(double));
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

double *solve_system(double **J, double *b, int n) {
    forward_elimination(J, b, n);
    return backward_substitution(J, b, n);
}

double* x_delta_by_gauss(double** J, double* b, int n) {
    solve_system(J,b,n);
}

void iter_newton(double* X, int n, double imax, double tol) {
    for (int i = 0; i < (int)imax; i++) {
        double** J = jacobian_exercise(X[0], X[1], X[2]);
        double* Y = function_exercise(X[0], X[1], X[2]);
        double* dX = x_delta_by_gauss(J, Y, n);
        for (int j = 0; j < n; j++) {
            X[j] -= dX[j];
        }
        double norm = 0;
        for (int j = 0; j < n; j++) {
            norm += dX[j] * dX[j];
        }
        norm = sqrt(norm);
        if (norm < tol) {
            printf("Converged.\n");
            break;
        }
        for (int j = 0; j < n; j++) {
            free(J[j]);
        }
        free(J);
        free(Y);
        free(dX);
    }
}

int main() {
    double X_0[3] = {100, 200, 3};
    iter_newton(X_0, 3, 10000, 0.00001);
    printf("Final solution: [%f, %f, %f]\n", X_0[0], X_0[1], X_0[2]);
    return 0;
}
