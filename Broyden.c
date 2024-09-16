#include <stdio.h>

double fabs(double x) {
    return (x < 0) ? -x : x;
}

double sqrt(double x) {
    if (x == 0 || x == 1) {
        return x;
    }

    double precision = 1e-6; 
    double guess = x / 2;

    while (fabs(guess * guess - x) > precision) {
        guess = (guess + x / guess) / 2;
    }

    return guess;
}

double fs(double x, double y) {
    return x + 2*y - 2;
}

double fs2(double x, double y) {
    return x*x + 4*y*y - 4;
}

void Js(double x, double y, double *J) {
    J[0] = 1;
    J[1] = 2;
    J[2] = 2;
    J[3] = 16;
}

int broyden_good(double *x, double *y, double tol, int maxIters) {
    int steps_taken = 0;
    double f[2], J[4];
    double s[2], z[2];
    double norm_f;

    f[0] = fs(*x, *y);
    f[1] = fs2(*x, *y);
    Js(*x, *y, J);

    while ((norm_f = sqrt(f[0]*f[0] + f[1]*f[1])) > tol && steps_taken < maxIters) {
        double det = J[0]*J[3] - J[1]*J[2];
        s[0] = (-J[3]*f[0] + J[1]*f[1]) / det;
        s[1] = (J[2]*f[0] - J[0]*f[1]) / det;

        *x += s[0];
        *y += s[1];

        double newf[2] = {fs(*x, *y), fs2(*x, *y)};
        z[0] = newf[0] - f[0];
        z[1] = newf[1] - f[1];

        J[0] += (z[0]*s[0] - J[0]*s[0]*s[0] - J[1]*s[0]*s[1]) / (s[0]*s[0] + s[1]*s[1]);
        J[1] += (z[0]*s[1] - J[0]*s[0]*s[1] - J[1]*s[1]*s[1]) / (s[0]*s[0] + s[1]*s[1]);
        J[2] += (z[1]*s[0] - J[2]*s[0]*s[0] - J[3]*s[0]*s[1]) / (s[0]*s[0] + s[1]*s[1]);
        J[3] += (z[1]*s[1] - J[2]*s[0]*s[1] - J[3]*s[1]*s[1]) / (s[0]*s[0] + s[1]*s[1]);

        f[0] = newf[0];
        f[1] = newf[1];
        steps_taken++;
    }

    return steps_taken;
}

int main() {
    double x0 = -1, y0 = -20;
    double tol = 1e-5;
    double maxIters = 100000;
    double n = broyden_good(&x0, &y0, tol, maxIters);
    printf("iterations: %lf\nx and y: %lf %lf\n", n, x0, y0);
    return 0;
}

