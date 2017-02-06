#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

void getHermiteWeights(double * x, double * w, int n);

double deltaFunction(int i, int j)
{
    if (i == j)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double solveGaussianHermiteQuadrature(int n,
                                      double (*Hx1)(double x1),
                                      double (*Hy1)(double y1),
                                      double (*Hx2)(double x2),
                                      double (*Hy2)(double y2),
                                      double (*Hx3)(double x3),
                                      double (*Hy3)(double y3),
                                      double (*Hx4)(double x4),
                                      double (*Hy4)(double y4),
                                      double (*V)(double x1 , double x2, double y1, double y2))
{
    /*
     * Solver for Gaussian Quadrature. Assumes that the integrand is a function(or class instance), that when called upon with
     * values x[i], x[j], x[k], x[l], will return the integrand value that is required.
     */

    double * x = new double[n]; // Integration points
    double * w = new double[n]; // Weights
    getHermiteWeights(x, w, n);
    double integralValue = 0;
    int i, j, k, l;

    for (i = 0; i < n; i++) // Loops over all x1
    {
        for (j = 0; j < n; j++) // Loops over all y1
        {
            for (k = 0; k < n; k++) // Loops over all x2
            {
                for (l = 0; l < n; l++) // Loops over all y2
                {
//                    CHANGE TO TAKE FOUR FUNCTIONS?
//                    integralValue += w[i]*w[j]*w[k]*w[l]*Integrand(x[i],x[j],x[k],x[l]);
                    integralValue += w[i]*w[j]*w[k]*w[l] * Hx1(x[i])*Hy1(x[j]) * Hx2(x[k])*Hy2(x[l]) * V(x[i],x[j],x[k],x[l]) * Hx3(x[i])*Hy3(x[j]) * Hx4(x[k])*Hy4(x[l]);
                }
            }
        }
    }

    delete [] x;
    delete [] w;

    return integralValue;
}

void getHermiteWeights(double * x, double * w, int n)
{
    /*
     * From MHJ lecture notes. Should probably check this one out closer.
     */
    int i, its, j, m;
    double p1,p2,p3,pp,z,z1;
    double Epsilon = 3.0e-14, PIM4 = 0.7511255444649425;
    int MaxIterations = 10;
    m=(n+1)/2;
    for (i=1;i<=m;i++)
    {
        if (i == 1)
        {
            z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
        }
        else if (i == 2)
        {
            z -= 1.14*pow((double)n,0.426)/z;
        }
        else if (i == 3)
        {
            z=1.86*z-0.86*x[0];
        }
        else if (i == 4)
        {
            z=1.91*z-0.91*x[1];
        }
        else
        {
            z=2.0*z-x[i-3];
        }
        for (its=1;its<=MaxIterations;its++)
        {
            p1=PIM4;
            p2=0.0;
            for (j=1;j<=n;j++)
            {
                p3=p2;
                p2=p1;
                p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
            }
            pp=sqrt((double)2*n)*p2;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= Epsilon) break;
        }
        if (its > MaxIterations) cout << "too many iterations in Hermite quadrature" << endl;
        x[i-1]=z;
        x[n-i] = -z;
        w[i-1]=2.0/(pp*pp);
        w[n-i]=w[i-1];
    }
}
