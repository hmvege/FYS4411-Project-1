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
//#include "omp.h"
#include "/usr/local/Cellar/llvm/3.9.1/lib/clang/3.9.1/include/omp.h" // To make omp.h usable on mac

#include "unittests.h"
//#include <time.h>

using namespace std;

double potentialV(double x1, double x2, double y1, double y2);
void printMatrix(double ** A, int N);
void setMatrixZero(double ** A, int N);

int main(int nargs, char *args[])
{
    int NElectrons  = 6; // Should be magic number: 2, 6, 12, 20
    int maxShell    = 6;


    int maxThreadNumber = omp_get_max_threads();
    int maxProcessorNumber = omp_get_num_procs();
    omp_set_num_threads(maxThreadNumber);
    cout << "  The number of processors available = " << maxProcessorNumber << endl ;
    cout << "  The number of threads available    = " << maxThreadNumber <<  endl;

    //    int intPoints   = atoi(args[1]);
//    int intPoints = 4;

    quantumDot QMDot(NElectrons, maxShell);
//    QMDot.setPotential(potentialV);
//    QMDot.setupInteractionMatrix(intPoints);
    QMDot.setupInteractionMatrixPolar();
//    QMDot.printInteractionMatrix(pow(maxShell,4));
    QMDot.runHartreeFock();

    /*
     * TODO:
     * - Add parallelization to integral-finder
     * - Make it so that hw is a function of omega, the basic energy and that omega is set from the main settup of program(should go into basis)
     * - Add possibility for looping over several electrons
     * - Add write-to-file capability
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
//            cout << A[i*N+j] << " "; // For contigious memory management
        }
        cout << endl;
    }
}

double potentialV(double x1, double x2, double y1, double y2)
{
    double eps = 1e-16;
//    double divisor = 2.0*((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    double divisor = 1.0*((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    if  (divisor > eps)
//        return 4.0/(sqrt(divisor));
        return 1.0/(sqrt(divisor));
    else
        return 0.0;
//    return 1.0/sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

void setMatrixZero(double ** A, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = 0;
        }
    }
}
