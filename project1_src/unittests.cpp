#include <iostream>
#include <iomanip>
#include <time.h>
#include "hermitepolynomials.h"
#include <cmath>
#include <vector>

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::vector;


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
//    int wrong = 0;
    double eps = 1e-15;
    for (int i = 0; i < maxHermite + 1; i++)
    {
        for (double x = 0; x < 10; x = x + 0.1)
        {
            cout << abs(Hermite.getPolynomRecursive(x,i) - Hermite.getPolynom(i)(x)) << endl;
            if (abs(Hermite.getPolynomRecursive(x,i) - Hermite.getPolynom(i)(x)) < eps)
            {
                correct++;
            }
            else
            {
                overflowX = x;
                overflowN = i;
//                wrong++;
                break;
            }
        }
    }
    cout << "Correct Hermite polynoms: " << correct << endl;
//    cout << "Wrong Hermite polynoms:   " << wrong << endl;
    if (overflowN > 0)
    {
        cout << "Most likely overflow encountered at Hermite polynomial N=" << overflowN << " and value x=" << overflowX << endl;
    }

    return 0;
}
