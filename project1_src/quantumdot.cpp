#include "quantumdot.h"
#include "hermitepolynomials.h"
#include "basis.h"
#include "functions.h"

quantumDot::quantumDot(int newNElectrons, int newMaxShell)
{
    N_Electrons = newNElectrons;
    maxShell = newMaxShell;
    basis.initializeBasis(maxShell);
    N_SPS = basis.getTotalParticleNumber();
}

void quantumDot::setPotential(double (*newV)(double x1, double x2, double y1, double y2))
{
    V = newV;
}

void quantumDot::setupInteractionMatrix()
{
    /*
     * Setting up the interaction matrix with contigious memory
     */
    interactionMatrix = new double[(int) pow(N_SPS,4)];

    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        double (*alphaHerm_x)(double) = hermite.getPolynom(basis.getState(alpha)->getN_x());
        double (*alphaHerm_y)(double) = hermite.getPolynom(basis.getState(alpha)->getN_y());
        double alphaConst = basis.getState(alpha)->normConst();
        for (int beta = 0; beta < N_SPS; beta++)
        {
            double (*betaHerm_x)(double) = hermite.getPolynom(basis.getState(beta)->getN_x());
            double (*betaHerm_y)(double) = hermite.getPolynom(basis.getState(beta)->getN_y());
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

                    interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)] = solveGaussianHermiteQuadrature(maxShell,
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
                    // state_alpha*alpha_norm * etc..
                    // May need to do some change of variable
                }
            }
        }
    }
}

void quantumDot::runHartreeFock()
{

}
