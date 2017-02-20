#include "hartreefock.h"
#include <armadillo>

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

int HartreeFock::runHF()
{
    // For loop will run till max HF iteration is reached, or we get a convergence for epsilon(NOT IMPLEMENTED YET)
    int maxHFiterations = 1000;
    int HFiteration = 0; // To compare with for later
    arma::vec oldEnergies = arma::zeros<arma::vec>(N_SPS);
    double lambda = 1e-8;

    std::cout << "Starting Hartree-Fock" << std::endl;
    for (HFiteration; HFiteration < maxHFiterations; HFiteration++)
    {
        // Setting up HFMatrix
        arma::mat HFMatrix = arma::zeros<arma::mat>(N_SPS,N_SPS);

        // Retrieving HF elements
        for (int alpha = 0; alpha < N_SPS; alpha++)
        {
            for (int beta = 0; beta < N_SPS; beta++)
            {
                double HFElement = 0;
                for (int gamma = 0; gamma < N_SPS; gamma++)
                {
                    for (int delta = 0; delta < N_SPS; delta++)
                    {
                        HFElement += densityMatrix(gamma,delta) * interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)];
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

        // Finding eigenvalues & eigenvectors
        arma::vec singleParticleEnergies;
        arma::eig_sym(singleParticleEnergies, C, HFMatrix);
        densityMatrix.zeros(); // Setting densityMatrix back to zeros only

        // Updating the density matrix
        updateDensityMatrix();

        // When if we have convergence, rather than checking everyone, find the max element
        if ((arma::sum(arma::abs(singleParticleEnergies - oldEnergies))/N_SPS) < lambda)
        {
            std::cout << "Done after " << HFiteration << " Hartree Fock iteration."<< std::endl;
            std::cout << singleParticleEnergies << std::endl;
            return 0;
        }
        else
        {
            oldEnergies = singleParticleEnergies;
        }

        if (HFiteration % 1000 == 0)
        {
            std::cout << HFiteration << std::endl;
        }
    }
    if (HFiteration == maxHFiterations)
    {
        std::cout << "Max HF iterations reached." << std::endl;
    }
    return -1; // Should probably remove this at some point...
}
