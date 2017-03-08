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
    int NElectronsArray[NElectronArrElems] = {2,6,12,20}; // Should be magic number: 2, 6, 12, 20. HAVE DONE: 2,6,PARTIALLY 12 & 20.NEED TO DO MORE ON 2 and 6 for higher shells
    int startShell          = 3;
    int maxShell            = 7;
    int maxHFIterations     = 200;
    double omega            = 1.0;
    double epsilon          = 1e-10;
    std::string filename    = "../output/HF_results";

    clock_t setupStart, setupFinish;
    setupStart = clock();

//    double HFEnergyPrev     = 0;
//    double HFEnergy         = 0;
//    double HFEnergyMaxDifference = 1e-4;

    for (int i = 1; i < NElectronArrElems-2; i++)
    {
        for (int shells = startShell; shells < maxShell; shells++)
        {
            quantumDot QMDot(NElectronsArray[i], shells, omega);
            QMDot.setupInteractionMatrixPolar();
            QMDot.setHFLambda(epsilon);
            QMDot.runHartreeFock(maxHFIterations);
            QMDot.storeResults(filename);
//            if (fabs(HFEnergy - HFEnergyPrev) < HFEnergyMaxDifference)
//            {
//                break;
//            }
//            else
//            {
//                HFEnergyPrev = HFEnergy;
//            }
        }
    }
    setupFinish = clock();
    cout << "Matrix setup complete. Setup time: " << ((setupFinish - setupStart)/((double)CLOCKS_PER_SEC)) << endl;
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
