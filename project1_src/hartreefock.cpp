#include "hartreefock.h"
#include <armadillo>
#include <iomanip>
#include <ctime>

using std::cout;
using std::endl;
using std::setprecision;

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
//    double lambda = 1e-5; // TODO: Add this as a user setting in QuantumDot class

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
                        // ADD M QM-NUMBER TEST HERE!
                        HFElement += densityMatrix(gamma,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
                    }
                }
                HFMatrix(alpha,beta) = HFElement;
//                if (alpha == beta) // Delta function
//                {
//                    HFMatrix(alpha,beta) += basis->getState(alpha)->getEnergy();
//                }
            }
            HFMatrix(alpha,alpha) += basis->getState(alpha)->getEnergy(); // Instead of an if-test inside the nested for-loop
        }
        loopFinish = clock();

        // Finding eigenvalues & eigenvectors
//        arma::vec singleParticleEnergies;
        eigStart = clock();
        singleParticleEnergies.zeros();
        arma::eig_sym(singleParticleEnergies, C, HFMatrix);
        densityMatrix.zeros(); // Setting densityMatrix back to zeros only
        eigFinish = clock();

        // Updating the density matrix
        updateDensityMatrix();

        minimaStart = clock();
        // When if we have convergence, rather than checking everyone, find the max element
        if ((arma::sum(arma::abs(singleParticleEnergies - oldEnergies)))/N_SPS < lambda)
        {
            SPS_Energies = singleParticleEnergies;
            writeToFile(singleParticleEnergies, C);
//            cout << singleParticleEnergies << endl;
            cout << "Done after " << HFiteration << " Hartree Fock iteration." << endl;
            break;
        }
        else
        {
            oldEnergies = singleParticleEnergies;
        }
        minimaFinish = clock();

        // Checking progress of program
        if (HFiteration % 1000 == 0)
        {
            cout << HFiteration << endl;
        }

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
        SPS_Energies = singleParticleEnergies;
        writeToFile(singleParticleEnergies, C);
//        cout << singleParticleEnergies << endl;
        cout << "Max HF iterations reached." << endl;
    }
    return 0; // Should probably remove this at some point...
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

    for (int i = 0; i < N_Electrons; i++)
    {
        energy += SPS_Energies(i);

                // Straight-forward method
//                energy += C(i,alpha)*C(i,beta)*basis->getState(alpha)->getEnergy();
//                for (int gamma = 0; gamma < N_SPS; gamma++)
//                {
//                    for (int delta = 0; delta < N_SPS; delta++)
//                    {
//                        for (int j = 0; j < N_Electrons; j++)
//                        {
//                            energy += 0.5*C(i,alpha)*C(j*beta)*C(i,gamma)*C(j,delta)*(interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
//                        }
//                    }
//                }
//                for (int gamma = 0; gamma < N_SPS; gamma++)
//                {
//                    for (int delta = 0; delta < N_SPS; delta)
//                    {
//                        energy += 0.5*densityMatrix(alpha, beta)*densityMatrix(gamma,delta)*(interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]);
//                    }
//                }
                // Simple summation from manipulating terms
//                energy += C(i,alpha)*C(i,beta)*SPS_Energies(alpha);

//            }
//        }
    }
    cout << "Energy: " << energy << endl;
}
