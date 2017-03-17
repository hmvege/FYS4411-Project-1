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
#include <mpi.h> // For mac
//#include "mpi/mpi.h" // For ubuntu

using namespace std;

int main(int numberOfArguments, char *cmdLineArguments[])
{
    int NElectronArrElems   = 4;
    int NElectronsArray[NElectronArrElems]; // Ugly setup
    NElectronsArray[0]      = 2;
    NElectronsArray[1]      = 6;
    NElectronsArray[2]      = 12;
    NElectronsArray[3]      = 20;

    int magicNumberIndex    = 3; // 0,1,2,3
    int startShell          = 4;
    int maxShell            = 10;
    int maxHFIterations     = 500;
    double omega            = 0.28;
    double epsilon          = 1e-10;
    std::string filename    = "../output2/HF_results";

    clock_t setupStart, setupFinish;
    setupStart = clock();

//    int numprocs, processRank;
//    MPI_Init (&numberOfArguments, &cmdLineArguments);
//    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
//    MPI_Comm_rank (MPI_COMM_WORLD, &processRank);

//    for (int i = magicNumberIndex; i < magicNumberIndex+1; i++)
    for (int i = 0; i < NElectronArrElems; i++)
    {
        for (int shells = startShell; shells < maxShell; shells++)
        {
            quantumDot QMDot(NElectronsArray[i], shells, omega);
            QMDot.setupInteractionMatrixPolar();
//            QMDot.setupInteractionMatrixPolarParalell(numprocs, processRank);
            QMDot.setHFLambda(epsilon);
//            QMDot.setOmega(0.5);
            QMDot.runHartreeFock(maxHFIterations);
            QMDot.storeResults(filename);
        }
    }

//    MPI_Finalize();

    setupFinish = clock();
    cout << "Program complete. Time used: " << ((setupFinish - setupStart)/((double)CLOCKS_PER_SEC)) << endl;
    /*
     * TODO:
     * [x] Clean up code quantumdot.cpp
     * [x] Clean up hartreefock.cpp
     * [x] Precalculate the antisymmetric integrals
     * [x] Add possibility for looping over several electrons(easy)
     * [ ] Add parallelization to integral-finder
     * [ ] Compare with unperturbed energy - unit test
     * [ ] Compare degeneracies before and afte
     * [ ] Add write-to-file capability(easy). For future project
     */
    return 0;
}
