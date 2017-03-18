#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "unittests.h"
#include "singlestate.h"
#include "basis.h"
#include "hermitepolynomials.h"
#include "functions.h"
#include "quantumdot.h"
#include "hartreefock.h"
#include <mpi.h> // For mac
//#include "mpi/mpi.h" // For ubuntu
//#include <libiomp/omp.h>

using namespace std;

int main(int numberOfArguments, char *cmdLineArguments[])
{
    testOrthogonality(numberOfArguments, cmdLineArguments);
    exit(1);
    cout<<"error"<<endl;
    int NElectronArrElems   = 4;
    int NElectronsArray[NElectronArrElems]; // Ugly setup
    NElectronsArray[0]      = 2;
    NElectronsArray[1]      = 6;
    NElectronsArray[2]      = 12;
    NElectronsArray[3]      = 20;

    int magicNumberIndex    = 3; // 0,1,2,3
    int startShell          = 3;
    int maxShell            = 6;
    int maxHFIterations     = 500;
    int noOmegas            = 3;
    double omega[noOmegas];
    omega[0]                = 1.0;
    omega[1]                = 0.5;
    omega[2]                = 0.28;
    double epsilon          = 1e-10;
    std::string filename    = "../output3/HF_results";

    clock_t setupStart, setupFinish;
    setupStart = clock();

    int numprocs, processRank;
    MPI_Init (&numberOfArguments, &cmdLineArguments);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &processRank);

    if (processRank == 0) { cout << "Starting up..." << endl; }

//    for (int i = magicNumberIndex; i < magicNumberIndex+1; i++)
    for (int i = 0; i < NElectronArrElems; i++)
    {
        for (int shells = startShell; shells < maxShell; shells++)
        {
            quantumDot QMDot(NElectronsArray[i], shells, omega[0]);
            if (false == checkElectronShellNumber(QMDot.getN_SPS(), QMDot.getN_Electrons())) { continue; }
            QMDot.initializeHF();
//            QMDot.setupInteractionMatrixPolar();
            QMDot.setupInteractionMatrixPolarParalell(numprocs, processRank);
//            QMDot.setupEmptyInteractionMatrix();
            if (processRank == 0)
            {
                for (int j = 1; j < noOmegas; j++)
                {
                    QMDot.setOmega(omega[j]);
                    QMDot.setHFLambda(epsilon);
                    QMDot.runHartreeFock(maxHFIterations);
                    QMDot.storeResults(filename);
                }
            }
        }
    }

    MPI_Finalize();

    setupFinish = clock();
    if (processRank == 0)
    {
        cout << "Program complete. Time used: " << ((setupFinish - setupStart)/((double)CLOCKS_PER_SEC)) << endl;
    }
    /*
     * TODO:
     * [x] Clean up code quantumdot.cpp
     * [x] Clean up hartreefock.cpp
     * [x] Precalculate the antisymmetric integrals
     * [x] Add possibility for looping over several electrons(easy)
     * [x] Add parallelization to integral-setup
     * [ ] Compare with unperturbed energy - unit test
     * [ ] Compare degeneracies before and after
     * [ ] Add write-to-file capability(easy). For future project
     */
    return 0;
}
