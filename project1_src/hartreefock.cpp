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
    initializeCMatrix();
    // Setting up the density matrix
    updateDensityMatrix();
}

void HartreeFock::initializeCMatrix()
{
    /*
     * Initializing the C-matrix.
     */
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
    /*
     * Function for updating the density matrix
     */
    for (int gamma = 0; gamma < N_SPS; gamma++)
    {
        for (int delta = 0; delta < N_SPS; delta++)
        {
            double sum = 0;
            for (int i = 0; i < N_Electrons; i++)
            {
                sum += C(gamma,i)*C(delta,i);
            }
            densityMatrix(gamma,delta) = sum;
        }
    }
}

void HartreeFock::setInteractionMatrix(double * newInteractionMatrix)
{
    /*
     * Function for setting the interaction matrix to be used in the algorithm.
     */
    interactionMatrix = newInteractionMatrix;
}

void HartreeFock::updateHFMatrix(arma::mat &HFMatrix)
{
    /*
     * Function for updating the Hartree-Fock matrix.
     */
    HFMatrix.zeros();
    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        int alpha_ml = basis->getState(alpha)->getM();
        int alphaSpin = basis->getState(alpha)->getSpin();
        for (int beta = 0; beta < N_SPS; beta++)
        {
            int beta_ml = basis->getState(beta)->getM();
            int betaSpin = basis->getState(beta)->getSpin();

            // Spin and M conservation test
            if ((alpha_ml != beta_ml) || (alphaSpin != betaSpin)) { continue; }

            HFMatrix(alpha,beta) = calculateInnerHFMatrixElement(alpha, alpha_ml, alphaSpin, beta, beta_ml, betaSpin);
        }
        // Instead of an if-test inside the nested for-loop
        HFMatrix(alpha,alpha) += basis->getState(alpha)->getEnergyPolar();
    }
}

double HartreeFock::calculateInnerHFMatrixElement(int alpha, int alpha_ml, int alphaSpin, int beta, int beta_ml, int betaSpin)
{
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
                HFElement += densityMatrix(gamma,delta) * interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)];
            }
        }
    }
    return HFElement;
}

int HartreeFock::runHF(int maxHFIterations)
{
    /*
     * Hartree-Fock algorithm. For loop will run till max HF iteration is reached, or we get a convergence in energies.
     */
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

    for (int HFIteration = 0; HFIteration < maxHFIterations; HFIteration++)
    {
        // Setting up HFMatrix
        loopStart = clock();
        updateHFMatrix(HFMatrix);
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
//            cout << "Done after " << HFIteration << " Hartree Fock iteration." << endl;
            HFCounter = HFIteration;
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

    // Checking if maximum HF iterations have been reached
    if (HFCounter == 0)
    {
        HFCounter = maxHFIterations;
        cout << "    Max HF iterations reached." << endl;
    }

    // Printing out average time per loop element
//    cout << "Average time per main loop:                 " << setprecision(8) << mainLoopTime / (double) HFCounter   << " seconds" << endl;
//    cout << "Average time for solving eigenvalueproblem: " << setprecision(8) << eigTime / (double) HFCounter        << " seconds" << endl;
//    cout << "Average time per finding minima:            " << setprecision(8) << minimaTime / (double) HFCounter     << " seconds" << endl;

    SPS_Energies = singleParticleEnergies;

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

void HartreeFock::printSPEnergies()
{
    /*
     * Small function for printing out the single particle energies with their belonging quantum numbers
     */
    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        cout << "SP-Energy = " << SPS_Energies(alpha) << " N = " << basis->getState(alpha)->getN() << " M = " << basis->getState(alpha)->getM() << endl;
    }
    cout << endl;
}

void HartreeFock::writeToFile()
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
        file << setprecision(8) << SPS_Energies(i) << " # ";
        for (int j = 0; j < N_SPS; j++)
        {
            file << setprecision(8) << C(i,j); // AM I FILLING IN EIGENVECTORS HERE?
        }
        file << endl;
    }
    file.close();
    cout << filename << " written" << endl;
}

void HartreeFock::getHFEnergy(double &HFEnergyResults, int &HFIterationsResults)
{
    /*
     * Function for retrieving the Hartree-Fock ground state energy
     */
    double energy = 0;

    for (int i = 0; i < N_Electrons; i++)
    {
        energy += SPS_Energies(i);
        for (int j = 0; j<N_Electrons; j++)
        {
            for (int alpha = 0; alpha < N_SPS; alpha++)
            {
                for (int beta = 0; beta < N_SPS; beta++)
                {
                    for (int gamma = 0; gamma < N_SPS; gamma++)
                    {
                        for (int delta= 0; delta < N_SPS; delta++)
                        {
                            energy += - 0.5 * C(alpha,i) * C(beta,i) * C(gamma,j) * C(delta,j) * interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)];
//                            energy += - 0.5 * C(alpha,i) * C(beta,i) * C(gamma,j) * C(delta,j) * sqrtOmega * interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)];
                        }
                    }
                }
            }
        }
    }
    if (processRank == 0)
    {
        printf("HF Iterations = %4d | Electrons = %2d | Shells = %2d | Energy = %3.6f \n", HFCounter, N_Electrons, basis->getMaxShell(), energy);
    }
    HFEnergyResults = energy;
    HFIterationsResults = HFCounter;
}
