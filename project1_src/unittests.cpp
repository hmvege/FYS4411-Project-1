#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <vector>
#include "hermitepolynomials.h"
#include "gaussianhermitequadrature.h"
#include "functions.h"

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::vector;

double PI = 3.14159265358979323846;

int runHermiteTimer(int N, double X)
{
    /* Benchmarking for the Hermitian polynom retriever.
     * N: number of iterations
     * X: value to insert into the Hermitian polynoms
     */

    HermitePolynomials herm;

    clock_t startHardCoded, finishHardCoded;
    startHardCoded = clock();

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 11; j++)
        {
            herm.getPolynom(j)(X);
        }
    }
    finishHardCoded = clock();

    clock_t startRecursive, finishRecursive;
    startRecursive = clock();

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 11; j++)
        {
            herm.getPolynomRecursive(X,j);
        }
    }
    finishRecursive = clock();

    cout << "Hardcoded results: " << setw(10) << setprecision(10) << (finishHardCoded - startHardCoded)/((double)CLOCKS_PER_SEC/1000) << endl;
    cout << "Recursive results: " << setw(10) << setprecision(10) << (finishRecursive - startRecursive)/((double)CLOCKS_PER_SEC/1000) << endl;

    return 0;
}

int checkCorrectHermitePolynoms(int maxHermite)
{
    /*
     * Checks that the different methods are equal
     */

    if (maxHermite > 10)
    {
        maxHermite = 10;
    }

    HermitePolynomials Hermite;

    double overflowX = 0;
    int overflowN = 0;

    int correct = 0;
    double eps = 1e-15;
    for (int i = 0; i < maxHermite + 1; i++)
    {
        for (double x = 0; x < 10; x = x + 0.1)
        {
            if (abs(Hermite.getPolynomRecursive(x,i) - Hermite.getPolynom(i)(x)) < eps)
            {
                correct++;
            }
            else
            {
                overflowX = x;
                overflowN = i;
                break;
            }
        }
    }
    cout << "Correct Hermite polynoms: " << correct << endl;
    if (overflowN > 0)
    {
        cout << "Most likely overflow encountered at Hermite polynomial N=" << overflowN << " and value x=" << overflowX << endl;
    }

    return 0;
}

double testFunction(double x)
{
    return 4*x*x;
}

double emptyFunction(double x)
{
    return 1;
}

double emptyPotentialFunction(double x1, double x2, double x3, double x4)
{
    return 1;
}

double testIntegrand(double x1, double x2, double x3, double x4)
{
    return testFunction(x1);
}

void testIntegratorClass(int NIntPoints, double eps=1e-10)
{
    /*
     * Small function for testing my integration class.
     */
    double analyticalSolution = 2*PI*PI;

    GaussianHermiteQuadrature integrator;
    integrator.setIntegrand(testIntegrand);
    double numericalSolution = integrator.solve(NIntPoints);

    printf("  Integration points: %d \n", NIntPoints);
    printf("  Numerical:  %.10f \n", numericalSolution);
    printf("  Analytical: %.10f \n", analyticalSolution);

    if (fabs(analyticalSolution - numericalSolution) < eps)
    {
        printf("  TEST PASSED \n");
    }
    else
    {
        printf("  TEST FAILED \n");
    }
}

void testIntegratorFunction(int NIntPoints, double eps=1e-10)
{
    /*
     * Small function for testing my integration function
     */

    double analyticalSolution = 2*PI*PI;

    // Not a pretty method
    double numericalSolution = solveGaussianHermiteQuadrature(NIntPoints, testFunction, emptyFunction, emptyFunction, emptyFunction, emptyFunction, emptyFunction, emptyFunction, emptyFunction, emptyPotentialFunction);

    printf("  Integration points: %d \n", NIntPoints);
    printf("  Numerical:  %.10f \n", numericalSolution);
    printf("  Analytical: %.10f \n", analyticalSolution);

    if (fabs(analyticalSolution - numericalSolution) < eps)
    {
        printf("  TEST PASSED \n");
    }
    else
    {
        printf("  TEST FAILED \n");
    }
}
