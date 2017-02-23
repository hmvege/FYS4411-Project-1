#include "quantumdot.h"
#include "hermitepolynomials.h"
#include "basis.h"
#include "functions.h"
#include "hartreefock.h"
#include "Coulomb_Functions.h"
#include <ctime>
//#include <omp.h>

using std::cout;
using std::endl;

quantumDot::quantumDot(int newNElectrons, int newMaxShell, double newOmega)
{
    N_Electrons = newNElectrons;
    maxShell = newMaxShell;
    omega = newOmega;
    basis.initializeBasis(maxShell, newOmega);
    N_SPS = basis.getTotalParticleNumber();
    interactionMatrixLength = (int) pow(N_SPS,4);
    HF.initializeHF(N_Electrons, N_SPS, &basis);

    //    basis.printBasisSize();
//    basis.printBasis(true); // For printing and checking basis
}

void quantumDot::setPotential(double (*newV)(double x1, double x2, double y1, double y2))
{
    V = newV;
}

void quantumDot::setupInteractionMatrix(int integrationPoints)
{
    /*
     * Setting up the interaction matrix with contigious memory for cartesian coordinates
     */

    clock_t setupStart, setupFinish;
    setupStart = clock();

    interactionMatrix = new double[interactionMatrixLength];

    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        double (*alphaHerm_x)(double) = hermite.getPolynom(basis.getState(alpha)->getN_x());
        double (*alphaHerm_y)(double) = hermite.getPolynom(basis.getState(alpha)->getN_y());
        double alphaConst = basis.getState(alpha)->normConst();
//        double constA1 = normalizationConstant(basis.getState(alpha)->getN_x());
//        double constA2 = normalizationConstant(basis.getState(alpha)->getN_y());
        for (int beta = 0; beta < N_SPS; beta++)
        {
            double (*betaHerm_x)(double) = hermite.getPolynom(basis.getState(beta)->getN_x());
            double (*betaHerm_y)(double) = hermite.getPolynom(basis.getState(beta)->getN_y());
//            double constA3 = normalizationConstant(basis.getState(beta)->getN_x());
//            double constA4 = normalizationConstant(basis.getState(beta)->getN_y());
            double betaConst = basis.getState(beta)->normConst();
            for (int gamma = 0; gamma < N_SPS; gamma++)
            {
                double (*gammaHerm_x)(double) = hermite.getPolynom(basis.getState(gamma)->getN_x());
                double (*gammaHerm_y)(double) = hermite.getPolynom(basis.getState(gamma)->getN_y());
                double gammaConst = basis.getState(gamma)->normConst();
                for (int delta = 0; delta < N_SPS; delta++)
                {
                    double (*deltaHerm_x)(double) = hermite.getPolynom(basis.getState(delta)->getN_x());
                    double (*deltaHerm_y)(double) = hermite.getPolynom(basis.getState(delta)->getN_y());
                    double deltaConst = basis.getState(delta)->normConst();

                    interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)] = solveGaussianHermiteQuadrature(integrationPoints,
                                                                                                                alphaHerm_x,
                                                                                                                alphaHerm_y,
                                                                                                                betaHerm_x,
                                                                                                                betaHerm_y,
                                                                                                                gammaHerm_x,
                                                                                                                gammaHerm_y,
                                                                                                                deltaHerm_x,
                                                                                                                deltaHerm_y,
                                                                                                                V);
                    interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)] *= alphaConst*betaConst*gammaConst*deltaConst;
//                    interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)] *= constA1*constA2*constA3*constA4;
                    // state_alpha*alpha_norm * etc..
                    // May need to do some change of variable..?
                }
            }
        }
    }
    setupFinish = clock();
    cout << "Matrix setup complete." << endl;
    cout << "Setup time: " << ((setupFinish - setupStart)/((double)CLOCKS_PER_SEC)) << endl;
}

void quantumDot::setupInteractionMatrixPolar()
{
    /*
     * Setting up the interaction matrix with contigious memory for polar coordinates
     */
    int nonEmptyStatesCounter = 0;
    interactionMatrix = new double[interactionMatrixLength];

    clock_t setupStart, setupFinish;
    setupStart = clock();

    double interactionValue = 0;
    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        int n1 = basis.getState(alpha)->getN();
        int ml1 = basis.getState(alpha)->getM();
        for (int gamma = alpha; gamma < N_SPS; gamma++)
        {
            int n2 = basis.getState(gamma)->getN();
            int ml2 = basis.getState(gamma)->getM();
            for (int beta = 0; beta < N_SPS; beta++)
            {
                int n3 = basis.getState(beta)->getN();
                int ml3 = basis.getState(beta)->getM();
                for (int delta = 0; delta < N_SPS; delta++) // CHANGE TO START AT BETA?
                {
                    int n4 = basis.getState(delta)->getN();
                    int ml4 = basis.getState(delta)->getM();

                    if ((basis.getState(alpha)->getSpin() + basis.getState(gamma)->getSpin()) !=
                            (basis.getState(beta)->getSpin() + basis.getState(delta)->getSpin()))
                    {
                        interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] = 0;
                        interactionMatrix[index(gamma, alpha, delta, beta, N_SPS)] = 0;
                    }
//                    else if (ml1 + ml2 != ml3 + ml4)
//                    {
//                        interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] = 0;
//                    }
                    else
                    {
                        nonEmptyStatesCounter++;
                        interactionValue = Coulomb_HO(basis.omega, n1, ml1, n2, ml2, n3, ml3, n4, ml4);
                        interactionMatrix[index(alpha, gamma, beta, delta, N_SPS)] = interactionValue;
                        interactionMatrix[index(gamma, alpha, delta, beta, N_SPS)] = interactionValue;
                    }
//                    interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)] = Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4);
                }
            }
        }
    }
    setupFinish = clock();
    cout << "Matrix setup complete. Number of non-empty states: " << nonEmptyStatesCounter << endl;
    cout << "Setup time: " << ((setupFinish - setupStart)/((double)CLOCKS_PER_SEC)) << endl;
}

void quantumDot::setupInteractionMatrixFromFile(const std::string& filename) // NOT FULLY IMPLEMENTED
{
    //intMatrixFromFile = HF.loadInteractionMatrix(filename);
    std::ifstream file(filename);
    //    std::string line;
    std::vector<double> inputVectorN;
    std::vector<double> inputVectorM;
    int n, m = 0;

    while (file >> n >> m)
    {
        inputVectorN.push_back(n - 1);
        inputVectorM.push_back(m - 1);
    }
}

void quantumDot::runHartreeFock(int maxHFIteration)
{
    HF.setInteractionMatrix(interactionMatrix);
    HF.runHF(maxHFIteration);
}

void quantumDot::printInteractionMatrix(int NPrintPoints)
{
    for (int i = 0; i < NPrintPoints; i++)
    {
        double intMatrix = interactionMatrix[i];
        if (intMatrix != 0.0)
        {
            printf("%.10f \n", interactionMatrix[i]);
        }
    }
}
