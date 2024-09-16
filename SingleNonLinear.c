#include <stdio.h>
#include <stdlib.h>

int n;

double fabs(double x) {
    return (x < 0) ? -x : x;
}

double power(int i, double x) {
    double temp = 1;
    for (int j = 0; j < i; j++) {
        temp *= x;
    }
    return temp;
}

double f(double x, double coef[]) {
    double sum = 0;
    for (int i = 0; i < n + 1; i++) {
        sum += coef[i] * power(i, x);
    }
    return sum;
}

double f_dash(double x, double der[]) {
    double sum = 0;
    for (int i = 0; i < n + 1; i++) {
        sum += der[i] * power(i, x);
    }
    return sum;
}

double incremental_search(double x[], double dx, double coef[]) {
    double f1, f2;
    x[1] = x[0] + dx;
    f1 = f(x[0], coef);
    f2 = f(x[1], coef);

    while (f1 * f2 > 0) {
        x[0] = x[1];
        x[1] = x[0] + dx;
        f1 = f2;
        f2 = f(x[1], coef);

    }

    if (f1 * f2 <= 0) {
        return x[0];
    }
    return -1;
}

double newton(double x[], double coef[], double der[], double Eo) {
    double dx;
    int i = 1;
    while (fabs(f(x[0], coef)) > Eo && i < 1000) { // Ensure iteration limit to avoid infinite loop
        dx = -(f(x[0], coef) / f_dash(x[0],der));
        x[1] = x[0] + dx;
        x[0] = x[1];
        i++;
    }
    return x[1];
}

double bisection(double x[], double coef[], double Eo) {
    double a = x[0];
    double b = x[1];
    double c;
    int i = 0;
    while (fabs(b - a) > Eo && i < 100000) { 
        c = (a + b) / 2;
        if (f(c, coef) == 0)
            return c;
        else if (f(a, coef) * f(c, coef) < 0)
            b = c;
        else
            a = c;
        i++;
    }
    return (a + b) / 2;
}

double regula_falsi(double coeff[], double x[], double Eo) {
    double curr;
    while (fabs(f(x[0], coeff) - f(x[1], coeff)) > Eo) {
        double first = f(x[0], coeff), second = f(x[1], coeff);
        curr = (second * x[0] - first * x[1]) / (second - first);
        if (f(curr, coeff) == 0 || fabs(curr - x[0]) < Eo) {
            break;
        }
        x[1] = x[0];
        x[0] = curr;
    }
    return curr;
}

int main(int argc, char const *argv[]) {
    printf("Enter the degree of polynomial: ");
    scanf("%d", &n);
    double *coef = malloc((n + 1) * sizeof(double));
    double *der = malloc((n + 1) * sizeof(double));
    if (coef == NULL || der == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    for (int i = 0; i < n + 1; i++) {
        printf("Enter coefficient of x^%d: ", i);
        scanf("%lf", &coef[i]);
        der[i] = i * coef[i];
    }
    double x[2], x_1[2];
    printf("Enter initial guess for Newton's method: ");
    scanf("%lf", &x[0]);
    double dx = 0.001;
    double Eo = 0.000001;
    double root1 = incremental_search(x, dx, coef);
    double root2 = newton(x, coef, der, Eo);
    printf("Enter initial guess for Regula-Falsi method: ");
    scanf("%lf %lf", &x_1[0], &x_1[1]);
    double root3 = regula_falsi(coef, x_1, Eo);
    double root4;
    printf("Enter initial interval for Bisection method: ");
    scanf("%lf %lf", &x[0], &x[1]);
    root4 = bisection(x, coef, Eo);

    if (root1 != -1) {
        printf("Root found by incremental search: %lf\n", root1);
    } else {
        printf("No root found by incremental search.\n");
    }
    printf("Root found by Newton's method: %lf\n", root2);
    printf("Root found by Regula-Falsi method: %lf\n", root3);
    printf("Root found by Bisection method: %lf\n", root4);
    free(coef);
    free(der);
    return 0;
}