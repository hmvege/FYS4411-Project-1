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
//#include <mpi.h> // For mac
//#include "mpi/mpi.h" // For ubuntu

using namespace std;

int main(int nargs, char *args[])
{
    int NElectronArrElems   = 4;
    int NElectronsArray[NElectronArrElems]; // Ugly setup
    NElectronsArray[0]      = 2;
    NElectronsArray[1]      = 6;
    NElectronsArray[2]      = 12;
    NElectronsArray[3]      = 20;

    int startShell          = 3;
    int maxShell            = 7;
    int maxHFIterations     = 500;
    double omega            = 1.0;
    double epsilon          = 1e-10;
    std::string filename    = "../output/HF_results";

    clock_t setupStart, setupFinish;
    setupStart = clock();

    for (int i = 1; i < NElectronArrElems-2; i++)
    {
        for (int shells = startShell; shells < maxShell; shells++)
        {
            quantumDot QMDot(NElectronsArray[i], shells, omega);
            QMDot.setupInteractionMatrixPolar();
            QMDot.setHFLambda(epsilon);
            QMDot.runHartreeFock(maxHFIterations);
            QMDot.storeResults(filename);
        }
    }
    setupFinish = clock();
    cout << "Program complete. Time used: " << ((setupFinish - setupStart)/((double)CLOCKS_PER_SEC)) << endl;
    /*
     * TODO:
     * [x] Clean up code quantumdot.cpp
     * [x] Clean up hartreefock.cpp
     * [x] Precalculate the antisymmetric integrals
     * [x] Add possibility for looping over several electrons(easy)
     * [ ] Add parallelization to integral-finder
     * [ ] Add parallelization to HF matrix setup
     * [ ] Add write-to-file capability(easy)
     * [ ] Compare with unperturbed energy - unit test
     */
    return 0;
}
