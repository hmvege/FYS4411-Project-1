#include "hartreefock.h"
#include <armadillo>
#include <iomanip>
#include <ctime>
#include "functions.h"
#include "Coulomb_Functions.h"

using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

HartreeFock::HartreeFock()
{

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

void HartreeFock::updateHFMatrix(arma::mat &HFMatrix)
{
    HFMatrix.zeros();
    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        int alpha_ml = basis->getState(alpha)->getM();
        int alphaSpin = basis->getState(alpha)->getSpin();
        for (int beta = 0; beta < N_SPS; beta++)
        {
            int beta_ml = basis->getState(beta)->getM();
            double betaSpin = basis->getState(beta)->getSpin();
            // Spin and M conservation test
            if ((alpha_ml != beta_ml) || (alphaSpin != betaSpin)) { continue; }
            double HFElement = 0;
            for (int gamma = 0; gamma < N_SPS; gamma++)
            {
                int gamma_ml = basis->getState(gamma)->getM();
                double gammaSpin = basis->getState(gamma)->getSpin();
                for (int delta = 0; delta < N_SPS; delta++)
                {
                    int delta_ml = basis->getState(delta)->getM();
                    double deltaSpin = basis->getState(delta)->getSpin();
                    // Total spin and angular momentum conservation test
                    if ((alpha_ml + gamma_ml == beta_ml + delta_ml) && (alphaSpin + gammaSpin == betaSpin + deltaSpin))
                    {
                        // Brute-forcey method, not removed the extra deltaFunction for spin
//                        HFElement += densityMatrix(gamma,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)]*deltaFunction(alphaSpin,betaSpin)*deltaFunction(gammaSpin,deltaSpin) - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]*deltaFunction(alphaSpin,deltaSpin)*deltaFunction(gammaSpin, betaSpin));
                        HFElement += densityMatrix(gamma,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)]*deltaFunction(gammaSpin,deltaSpin) - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]*deltaFunction(alphaSpin,deltaSpin)*deltaFunction(gammaSpin, betaSpin));
                    }
                }
            }
            HFMatrix(alpha,beta) = HFElement;
        }
        // Instead of an if-test inside the nested for-loop
        HFMatrix(alpha,alpha) += basis->getState(alpha)->getEnergyPolar();
    }
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
        updateHFMatrix(HFMatrix);

        //        HFMatrix.zeros();
//        for (int alpha = 0; alpha < N_SPS; alpha++)
//        {
//            int alpha_ml = basis->getState(alpha)->getM();
//            int alphaSpin = basis->getState(alpha)->getSpin();

//            for (int beta = 0; beta < N_SPS; beta++)
//            {
//                int beta_ml = basis->getState(beta)->getM();
//                double betaSpin = basis->getState(beta)->getSpin();

//                // Spin and M conservation test
//                if ((alpha_ml != beta_ml) || (alphaSpin != betaSpin))
//                {
//                    continue;
//                }

//                double HFElement = 0;
//                for (int gamma = 0; gamma < N_SPS; gamma++)
//                {
//                    int gamma_ml = basis->getState(gamma)->getM();
//                    double gammaSpin = basis->getState(gamma)->getSpin();

//                    for (int delta = 0; delta < N_SPS; delta++)
//                    {
//                        int delta_ml = basis->getState(delta)->getM();
//                        double deltaSpin = basis->getState(delta)->getSpin();

//                        // Total spin and angular momentum conservation test
//                        if ((alpha_ml + gamma_ml == beta_ml + delta_ml) && (alphaSpin + gammaSpin == betaSpin + deltaSpin))
//                        {
////                            HFElement += densityMatrix(gamma,delta) * (Coulomb_HO(basis->omega, alpha_n, alpha_ml, gamma_n, gamma_ml, beta_n, beta_ml, delta_n, delta_ml)*deltaFunction(alphaSpin,betaSpin)*deltaFunction(gammaSpin,deltaSpin) - Coulomb_HO(basis->omega, alpha_n, alpha_ml, gamma_n, gamma_ml, delta_n, delta_ml, beta_n, beta_ml)*deltaFunction(alphaSpin,deltaSpin)*deltaFunction(gammaSpin, betaSpin));
//                            HFElement += densityMatrix(gamma,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)]*deltaFunction(alphaSpin,betaSpin)*deltaFunction(gammaSpin,deltaSpin) - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]*deltaFunction(alphaSpin,deltaSpin)*deltaFunction(gammaSpin, betaSpin));
//                        }

//                    }
//                }
//                HFMatrix(alpha,beta) = HFElement;
//            }
//            // Instead of an if-test inside the nested for-loop
//            HFMatrix(alpha,alpha) += basis->getState(alpha)->getEnergyPolar();
//        }
        loopFinish = clock();

        // Finding eigenvalues & eigenvectors
        eigStart = clock();
        singleParticleEnergies.zeros();
        arma::eig_sym(singleParticleEnergies, C, HFMatrix); // dc faster for larger matrices, std default
        densityMatrix.zeros(); // Setting densityMatrix back to zeros only
        eigFinish = clock();

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
    cout<<SPS_Energies<<endl;
    writeToFile(singleParticleEnergies, C);

    return 0;
}

void HartreeFock::printHFMatrix(arma::mat HFMatrix)
{
    /*
     *  Small function for printing out the Hartree-Fock matrix
     */
    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        for (int beta = 0; beta < N_SPS; beta++)
        {
            double val = HFMatrix(alpha,beta);
            if (fabs(val) < (1e-15))
            {
                cout << setw(12) << 0;
            }
            else
            {
                cout << setw(12) << HFMatrix(alpha,beta);
            }
        }
    cout << endl;
    }
}

void HartreeFock::printSPEnergies(arma::vec singleParticleEnergies)
{
    /*
     * Small function for printing out the single particle energies with their belonging quantum numbers
     */
    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        cout << "SP-Energy = " << singleParticleEnergies(alpha) << " N = " << basis->getState(alpha)->getN() << " M = " << basis->getState(alpha)->getM() << endl;
    }
    cout << endl;
}

void HartreeFock::writeToFile(arma::vec eigVals, arma::mat eigVecs)
{
    /*
     * Function for writing eigenvectors and the egeinvalue to file
     * Setup:
     * eigenValue # eigenvector-elements
     */
    std::ofstream file;
    std::string filename = "output_NSPS" + std::to_string(N_SPS) + "_Electrons" + std::to_string(N_Electrons) + ".dat";
    file.open(filename);
    for (int i = 0; i < N_SPS; i++)
    {
        file << setprecision(8) << eigVals(i) << " # ";
        for (int j = 0; j < N_SPS; j++)
        {
            file << setprecision(8) << eigVecs(i,j); // AM I FILLING IN EIGENVECTORS HERE?
        }
        file << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}

void HartreeFock::getHFEnergy()
{
    /*
     * Function for retrieving the Hartree-Fock ground state energy
     */

    double energy = 0;
    double alphaSpin = 0;
    double betaSpin = 0;
    double gammaSpin = 0;
    double deltaSpin = 0;

//    // Test for checking unity of C matrix
//    if (fabs(arma::sum(arma::sum(C*arma::trans(C))) - N_SPS) < lambda)
//    {
//        cout << "tr(C*C^T) = " << arma::trace(C*arma::trans(C)) << endl;
//        cout << "C is unitary" << endl;
//    }

    for (int i = 0; i < N_Electrons; i++)
    {
        energy += SPS_Energies(i);
        for (int j = 0; j<N_Electrons; j++)
        {
        for (int alpha = 0; alpha < N_SPS; alpha++)
        {
            alphaSpin = basis->getState(alpha)->getSpin();
            for (int beta = 0; beta < N_SPS; beta++)
            {
                betaSpin = basis->getState(beta)->getSpin();
                for (int gamma = 0; gamma < N_SPS; gamma++)
                {
                    gammaSpin = basis->getState(gamma)->getSpin();
                    for (int delta= 0; delta < N_SPS; delta++)
                    {
                        deltaSpin = basis->getState(delta)->getSpin();
                        energy += - 0.5 * C(alpha,i) * C(beta,i) * C(gamma,j) * C(delta,j) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)]*deltaFunction(alphaSpin,betaSpin)*deltaFunction(gammaSpin,deltaSpin) - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]*deltaFunction(alphaSpin,deltaSpin)*deltaFunction(gammaSpin,betaSpin));

                        //                        energy += - 0.5 * densityMatrix(alpha,gamma) * densityMatrix(beta,delta) * (interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)]*deltaFunction(alphaSpin,betaSpin)*deltaFunction(gammaSpin,deltaSpin) - interactionMatrix[index(alpha, gamma, delta, beta, N_SPS)]*deltaFunction(alphaSpin,deltaSpin)*deltaFunction(gammaSpin,betaSpin));
                        }
                    }
                }
            }
        }

    }

    cout << "Energy: " << energy << endl;
}
