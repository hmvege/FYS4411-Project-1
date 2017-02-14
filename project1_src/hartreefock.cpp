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
    for (int i = 0; i < N_SPS; i++)
    {
        C(i,i) = 1;
    }

    // Setting up the density matrix
    for (int gamma = 0; gamma < N_SPS; gamma++)
    {
        for (int delta = 0; delta < N_SPS; delta++)
        {
            double sum = 0;
            for (int i = 0; i < N_Electrons; i++)
            {
                sum += C(i,gamma)*C(i,delta);
            }
            densityMatrix(gamma,delta) = sum;
        }
    }
}

void HartreeFock::initializeHF(int NElectrons, int singleParticleStates)
{
    N_Electrons = NElectrons;
    N_SPS = singleParticleStates;

    densityMatrix = arma::zeros<arma::mat>(N_SPS,N_SPS);
    C = arma::zeros<arma::mat>(N_SPS,N_SPS);

    // Initializing C as a diagonal matrix(other initializations exists)
    for (int i = 0; i < N_SPS; i++)
    {
        C(i,i) = 1;
    }

    // Setting up the density matrix
    for (int gamma = 0; gamma < N_SPS; gamma++)
    {
        for (int delta = 0; delta < N_SPS; delta++)
        {
            double sum = 0;
            for (int i = 0; i < N_Electrons; i++)
            {
                sum += C(i,gamma)*C(i,delta);
            }
            densityMatrix(gamma,delta) = sum;
        }
    }
}


//HartreeFock::HartreeFock(int NElectrons, int singleParticleStates)
//{
//    N_Electrons = NElectrons;
//    N_SPS = singleParticleStates;

////    densityMatrix = new double * [N_SPS];
////    C = new double * [N_SPS];

//    densityMatrix = arma::zeros<arma::mat>(N_SPS,N_SPS);
//    C = arma::zeros<arma::mat>(N_SPS,N_SPS);

////    for (int i = 0; i < N_SPS; i++)
////    {
////        C[i] = new double[N_SPS];
////        densityMatrix[i] = new double[N_SPS];
////    }

//    // Initializing C as a diagonal matrix(other initializations exists)
//    for (int i = 0; i < N_SPS; i++)
//    {
////        C[i][i] = 1;
//        C(i,i) = 1;
//    }

//    // Setting up the density matrix
//    for (int gamma = 0; gamma < N_SPS; gamma++)
//    {
//        for (int delta = 0; delta < N_SPS; delta++)
//        {
//            double sum = 0;
//            for (int i = 0; i < N_Electrons; i++)
//            {
////                sum += C[i][gamma]*C[i][delta];
//                sum += C(i,gamma)*C(i,delta);
//            }
////            densityMatrix[gamma][delta] = sum;
//            densityMatrix(gamma,delta) = sum;
//        }
//    }
//}

HartreeFock::~HartreeFock()
{
//    for (int i = 0; i < N_Electrons; i++)
//    {
//        delete [] densityMatrix[i];
//        delete [] C[i];
//    }
//    delete [] densityMatrix;
//    delete [] C;
    delete [] interactionMatrix;
}

void HartreeFock::setInteractionMatrix(double * newInteractionMatrix)
{
    interactionMatrix = newInteractionMatrix;
}

int HartreeFock::runHF()
{
    // For loop will run till max HF iteration is reached, or we get a convergence for epsilon(NOT IMPLEMENTED YET)
    int maxHFiterations = 100;
    arma::vec oldEnergies = arma::zeros<arma::vec>(N_SPS);
    double lambda = 1e-8;

    for (int HFiteration = 0; HFiteration < maxHFiterations; HFiteration++)
    {
        // Setting up HFMatrix
        arma::mat HFMatrix = arma::zeros<arma::mat>(N_SPS,N_SPS);

//        double ** HFMatrix = new double * [N_SPS];
//        for (int i = 0; i < N_SPS; i++) {
//            HFMatrix[i] = new double[N_SPS];
//        }
//        for (int i = 0; i < N_SPS; i++) {
//            for (int j = 0; j < N_SPS; j++) {
//                HFMatrix[i][j] = 0;
//            }
//        }

        // retrieving HF elements
        for (int alpha = 0; alpha < N_SPS; alpha++)
        {
            for (int beta = 0; beta < N_SPS; beta++)
            {
                double HFElement = 0;
                for (int gamma = 0; gamma < N_SPS; gamma++)
                {
                    for (int delta = 0; delta < N_SPS; delta++)
                    {
//                        HFElement += C[gamma][delta] * interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)];
                        HFElement += C(gamma,delta) * interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)];
                    }
                }
                HFMatrix(alpha,beta) = HFElement;
//                HFMatrix[alpha][beta] = HFElement;
                if (alpha == beta) // Delta function
                {
//                    HFMatrix[alpha][beta] += basis.getState(alpha)->getEnergy();
                    HFMatrix(alpha,beta) += basis.getState(alpha)->getEnergy();
                }
            }
        }

        // Finding eigenvalues & eigenvectors
        arma::vec singleParticleEnergies;
        arma::eig_sym(singleParticleEnergies, C, HFMatrix);

        densityMatrix.zeros(); // Setting densityMatrix back to zeros only

        // Rerunning the
        for (int gamma = 0; gamma < N_SPS; gamma++)
        {
            for (int delta = 0; delta < N_SPS; delta++)
            {
                double sum = 0;
                for (int i = 0; i < N_Electrons; i++)
                {
                    sum += C(gamma,i) * C(delta,i);
                }
                densityMatrix(gamma,delta) = sum;
            }
        }

        if ((fabs(arma::sum(singleParticleEnergies) - arma::sum(oldEnergies))) < lambda)
        {
            std::cout << "Done" << std::endl;
            return 0;
        }
        else
        {
            oldEnergies = singleParticleEnergies;
        }

        // When if we have convergence, rather than checking everyone, find the max element

        // De-allocating memory
//        for (int i = 0; i < N_SPS; i++) {
//            delete [] HFMatrix[i];
//        }
//        delete [] HFMatrix;
    }
    return -1;
}
