#ifndef FUNCTION_H
#define FUNCTION_H

double deltaFunction(int i, int j);
int factorial(int n);
int index(int i, int j, int k, int l, int N);
double solveGaussianHermiteQuadrature(int n,
                                      double (*Hx1)(double x1),
                                      double (*Hy1)(double y1),
                                      double (*Hx2)(double x2),
                                      double (*Hy2)(double y2),
                                      double (*Hx3)(double x3),
                                      double (*Hy3)(double y3),
                                      double (*Hx4)(double x4),
                                      double (*Hy4)(double y4),
                                      double (*V)(double x1, double x2, double y1, double y2));
#endif // FUNCTION_H
