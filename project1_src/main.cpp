#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "singlestate.h"
#include "basis.h"
#include "hermitepolynomials.h"
#include "functions.h"
#include "quantumdot.h"
#include "hartreefock.h"
//#include <mpi.h>

using namespace std;

double potentialV(double x1, double x2, double y1, double y2);
void printMatrix(double ** A, int N);
void setMatrixZero(double ** A, int N);

int main(int nargs, char *args[])
{
    int NElectrons  = 6; // Should be magic number: 2, 6, 12, 20
    int maxShell    = 3;
    int maxHFIterations = 200;
    double omega = 1;
    double epsilon = 1e-10;

//    int maxThreadNumber = omp_get_max_threads();
//    int maxProcessorNumber = omp_get_num_procs();
//    omp_set_num_threads(maxThreadNumber);
//    cout << "  The number of processors available = " << maxProcessorNumber << endl ;
//    cout << "  The number of threads available    = " << maxThreadNumber <<  endl;

    quantumDot QMDot(NElectrons, maxShell, omega);

    QMDot.setupInteractionMatrixPolar();
    QMDot.setHFLambda(epsilon);
    QMDot.runHartreeFock(maxHFIterations);

    /*
     * TODO:
     * [x] Clean up code quantumdot.cpp
     * [x] Clean up hartreefock.cpp
     * [x] Precalculate the antisymmetric integrals
     * [ ] Add possibility for looping over several electrons(easy)
     * [ ] Add parallelization to integral-finder
     * [ ] Add parallelization to HF matrix setup
     * [ ] Add write-to-file capability(easy)
     */

    return 0;
}

void printMatrix(double ** A, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

double potentialV(double x1, double x2, double y1, double y2)
{
    double eps = 1e-16;
    double divisor = 1.0*((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    if  (divisor > eps)
        return 1.0/(sqrt(divisor));
    else
        return 0.0;
//    return 1.0/sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}
