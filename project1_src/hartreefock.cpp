#include "hartreefock.h"
#include <armadillo>
#include <iomanip>
#include <ctime>

using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

HartreeFock::HartreeFock()
{

}

HartreeFock::HartreeFock(int NElectrons, int singleParticleStates)
{
    N_Electrons = NElectrons;
    N_SPS = singleParticleStates;

    densityMatrix = arma::zeros<arma::mat>(N_SPS,N_SPS);
    C = arma::zeros<arma::mat>(N_SPS,N_SPS);

    // Initializing C as a diagonal matrix(other initializations exists)
    setCMatrix();

    // Setting up the density matrix
    updateDensityMatrix();
}

void HartreeFock::initializeHF(int NElectrons, int singleParticleStates, Basis *newBasis)
{
    N_Electrons = NElectrons;
    N_SPS = singleParticleStates;
    basis = newBasis;
    densityMatrix = arma::zeros<arma::mat>(N_SPS,N_SPS);
    C = arma::zeros<arma::mat>(N_SPS,N_SPS);

    // Initializing C as a diagonal matrix(other initializations exists)
    setCMatrix();

    // Setting up the density matrix
    updateDensityMatrix();
}

void HartreeFock::setCMatrix()
{
    for (int i = 0; i < N_SPS; i++)
    {
        C(i,i) = 1;
    }
}

HartreeFock::~HartreeFock()
{
    delete [] interactionMatrix;
}

void HartreeFock::updateDensityMatrix()
{
    for (int gamma = 0; gamma < N_SPS; gamma++)
    {
        for (int delta = 0; delta < N_SPS; delta++)
        {
            double sum = 0;
            for (int i = 0; i < N_Electrons; i++)
            {
                sum += C(gamma,i)*C(delta,i); // Is this right??
            }
            densityMatrix(gamma,delta) = sum;
        }
    }
}

void HartreeFock::setInteractionMatrix(double * newInteractionMatrix)
{
    interactionMatrix = newInteractionMatrix;
}

int HartreeFock::runHF(int maxHFIterations)
{
    // For loop will run till max HF iteration is reached, or we get a convergence for epsilon(NOT IMPLEMENTED YET)
    arma::vec oldEnergies = arma::zeros<arma::vec>(N_SPS);
    arma::vec singleParticleEnergies = arma::zeros<arma::vec>(N_SPS);
    arma::mat HFMatrix = arma::zeros<arma::mat>(N_SPS,N_SPS);

    // For timing functions
    double mainLoopTime = 0;
    double eigTime = 0;
    double minimaTime = 0;
    clock_t loopStart, loopFinish;
    clock_t eigStart, eigFinish;
    clock_t minimaStart, minimaFinish;

    cout << "Starting Hartree-Fock" << endl;
    int HFiteration = 0; // To compare with for later(I know this is risky business)
    for (HFiteration; HFiteration < maxHFIterations; HFiteration++)
    {
        cout << "HFiteration: " << HFiteration << endl; // Checking progress of program, usefull for short runs

        // Setting up HFMatrix
        loopStart = clock();
        HFMatrix.zeros();
        // Retrieving HF elements; TODO: PARALELLIZE HERE
        for (int alpha = 0; alpha < N_SPS; alpha++)
        {
            for (int beta = 0; beta < N_SPS; beta++)
            {
                // TODO: Add quantum-number conservation tests here! SPIN TEST HERE?!

                double HFElement = 0;
                for (int gamma = 0; gamma < N_SPS; gamma++)
                {
                    for (int delta = 0; delta < N_SPS; delta++)
                    {
//                        // ADD M QM-NUMBER TEST HERE!
//                        if (basis->getState(alpha)->getM() + basis->getState(gamma)->getM() == basis->getState(beta)->getM() + basis->getState(delta)->getM())
//                        {
//                            HFElement += densityMatrix(gamma,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
//                        }
                        HFElement += densityMatrix(gamma,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
                    }
                }
                HFMatrix(alpha,beta) = HFElement;
//                if (alpha == beta) // Delta function
//                {
//                    HFMatrix(alpha,beta) += basis->getState(alpha)->getEnergy();
//                }
            }
            HFMatrix(alpha,alpha) += basis->getState(alpha)->getEnergyPolar(); // Instead of an if-test inside the nested for-loop
        }
        loopFinish = clock();

//        // Printing HFMatrix
//        for (int alpha = 0; alpha < N_SPS; alpha++)
//        {
//            for (int beta = 0; beta < N_SPS; beta++)
//            {
//                double val = HFMatrix(alpha,beta);
//                if (fabs(val) < (1e-15))
//                {
//                    cout << setw(12) << 0;
//                }
//                else
//                {
//                    cout << setw(12) << HFMatrix(alpha,beta);
//                }
//            }
//            cout << endl;
//        }

        // Finding eigenvalues & eigenvectors
        eigStart = clock();
        singleParticleEnergies.zeros();
        arma::eig_sym(singleParticleEnergies, C, HFMatrix); // dc faster for larger matrices, std default
        densityMatrix.zeros(); // Setting densityMatrix back to zeros only
        eigFinish = clock();

//        for (int alpha = 0; alpha < N_SPS; alpha++)
//        {
//            cout << singleParticleEnergies(alpha) << " N=" << basis->getState(alpha)->getN() << " M="<< basis->getState(alpha)->getM() << " SPEnergy: " << basis->getState(alpha)->getEnergyPolar() << endl;
//        }
//        cout << endl;


        // Updating the density matrix
        updateDensityMatrix();

        // When if we have convergence, rather than checking everyone, find the max element
        minimaStart = clock();
        if ((arma::sum(arma::abs(singleParticleEnergies - oldEnergies)))/N_SPS < lambda)
        {
//            cout << singleParticleEnergies << endl;
            cout << "Done after " << HFiteration << " Hartree Fock iteration." << endl;
            break;
        }
        else
        {
            oldEnergies = singleParticleEnergies;
        }
        minimaFinish = clock();

        // Checking progress of program, usefull for long runs
//        if (HFiteration % 1000 == 0)
//        {
//            cout << HFiteration << endl;
//        }

        // Summing up time spent
        mainLoopTime += ((loopFinish - loopStart)/((double)CLOCKS_PER_SEC));
        eigTime += ((eigFinish - eigStart)/((double)CLOCKS_PER_SEC));
        minimaTime += ((minimaFinish - minimaStart)/((double)CLOCKS_PER_SEC));
    }

    // Printing out average time per loop element
    cout << "Average time per main loop:                 " << setprecision(8) << mainLoopTime / (double) HFiteration   << " seconds" << endl;
    cout << "Average time for solving eigenvalueproblem: " << setprecision(8) << eigTime / (double) HFiteration        << " seconds" << endl;
    cout << "Average time per finding minima:            " << setprecision(8) << minimaTime / (double) HFiteration     << " seconds" << endl;

    if (HFiteration == maxHFIterations)
    {
        cout << "Max HF iterations reached." << endl;
    }
    SPS_Energies = singleParticleEnergies;
    writeToFile(singleParticleEnergies, C);

    return 0;
}

void HartreeFock::writeToFile(arma::vec eigVals, arma::mat eigVecs)
{
    std::ofstream file;
    std::string filename = "output_NSPS" + std::to_string(N_SPS) + "_Electrons" + std::to_string(N_Electrons) + ".dat";
    file.open(filename);
    for (int i = 0; i < N_SPS; i++)
    {
        file << setprecision(8) << eigVals(i) << " # ";
//        for (int j = 0; j < N_SPS; j++)
//        {
//            file << setprecision(8) << eigVecs(i,j); // AM I FILLING IN EIGENVECTORS HERE?
//        }
        file << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}

void HartreeFock::getEnergies()
{
    double energy = 0;
    double matrixEnergies = 0;

    double alphaSpin = 0;
    double betaSpin = 0;
    double gammaSpin = 0;
    double deltaSpin = 0;

    // Test for checking unity of C matrix
    if (fabs(arma::sum(arma::sum(C*arma::trans(C))) - N_SPS) < lambda)
    {
        cout<< "C*C^T = " << arma::sum(arma::sum(C*arma::trans(C)))<<endl;
        cout << "C is unitary" << endl;
    }

    for (int i = 0; i < N_Electrons; i++)
    {
        for (int j = 0; j < N_Electrons; j++)
        {
            for (int alpha = 0; alpha < N_SPS; alpha++)
            {
                for (int beta = 0; beta < N_SPS; beta++)
                {
                    for (int gamma = 0; gamma < N_SPS; gamma++)
                    {
                        for (int delta= 0; delta < N_SPS; delta++)
                        {
                            energy += - 0.5 * C(alpha,i) * C(beta,j) * C(gamma,i) * C(delta,j) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
                        }
                    }
                }
            }
        }
    }

//    for (int i = 0; i < N_Electrons; i++)
//    {
//        energy += SPS_Energies(i);
//        for (int alpha = 0; alpha < N_SPS; alpha++)
//        {
//            alphaSpin = basis->getState(alpha)->getSpin();
//            for (int beta = 0; beta < N_SPS; beta++)
//            {
//                betaSpin = basis->getState(beta)->getSpin();
//                for (int gamma = 0; gamma < N_SPS; gamma++)
//                {
//                    gammaSpin = basis->getState(gamma)->getSpin();
//                    for (int delta= 0; delta < N_SPS; delta++)
//                    {
//                        deltaSpin = basis->getState(delta)->getSpin();
//                        if (((alphaSpin + gammaSpin) == (betaSpin + deltaSpin)) && (basis->getState(alpha)->getM() + basis->getState(gamma)->getM() == basis->getState(beta)->getM() + basis->getState(delta)->getM()))
////                        if (((alphaSpin + gammaSpin) == (betaSpin + deltaSpin)) && ((basis->getState(alpha)->getM() + basis->getState(beta)->getM()) == basis->getState(gamma)->getM() + basis->getState(delta)->getM()))
//                        {
//                            double en = - 0.5 * densityMatrix(alpha,gamma) * densityMatrix(beta,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
////                            if (fabs(en) > lambda)
////                            {
////                                cout << en << endl;
////                            }
//                            matrixEnergies += en;
////                            energy -= 0.5 * densityMatrix(alpha,gamma) * densityMatrix(beta,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
////                            energy += - 0.5 * C(alpha,i) * C(beta,j) * C(gamma,i) * C(delta,j) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
//                        }
//                    }
//                }
//            }
//        }
//    }

    cout << "Matrix Energies: " << matrixEnergies << endl;
    cout << "Energy: " << energy << endl;
}
