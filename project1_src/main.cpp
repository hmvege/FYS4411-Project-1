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

int main(int nargs, char *args[])
{
    int NElectronArrElems   = 2;
    int NElectronsArray[NElectronArrElems] = {12,20}; // Should be magic number: 2, 6, 12, 20. HAVE DONE: 2,6
    int startShell          = 4;
    int maxShell            = 12;
    int maxHFIterations     = 200;
    double omega            = 1.0;
    double epsilon          = 1e-10;
    std::string filename    = "../output/HF_results.txt";

//    int maxThreadNumber = omp_get_max_threads();
//    int maxProcessorNumber = omp_get_num_procs();
//    omp_set_num_threads(maxThreadNumber);
//    cout << "  The number of processors available = " << maxProcessorNumber << endl ;
//    cout << "  The number of threads available    = " << maxThreadNumber <<  endl;

    double HFEnergyPrev     = 0;
    double HFEnergy         = 0;
    double HFEnergyMaxDifference = 1e-4;

    for (int i = 0; i < NElectronArrElems; i++)
    {
        for (int shells = startShell; shells < maxShell; shells++)
        {
            quantumDot QMDot(NElectronsArray[i], shells, omega);
            QMDot.setupInteractionMatrixPolar();
            QMDot.setHFLambda(epsilon);
            HFEnergy = QMDot.runHartreeFock(maxHFIterations);
            QMDot.storeResults(filename);
            if (fabs(HFEnergy - HFEnergyPrev) < HFEnergyMaxDifference)
            {
                break;
            }
            else
            {
                HFEnergyPrev = HFEnergy;
            }
        }
    }

    /*
     * TODO:
     * [x] Clean up code quantumdot.cpp
     * [x] Clean up hartreefock.cpp
     * [x] Precalculate the antisymmetric integrals
     * [x] Add possibility for looping over several electrons(easy)
     * [ ] Add parallelization to integral-finder
     * [ ] Add parallelization to HF matrix setup
     * [ ] Add write-to-file capability(easy)
     * [ ] Compare with unperturbed energy - unittest
     */

    return 0;
}
