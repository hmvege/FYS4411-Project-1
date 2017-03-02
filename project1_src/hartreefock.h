#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <armadillo>
#include "basis.h"
#include "hermitepolynomials.h"
#include "Coulomb_Functions.h"

class HartreeFock
{
private:
    int N_SPS;
    int N_Electrons;
    double lambda = 1e-5;
    double * interactionMatrix;
    arma::vec SPS_Energies;
    arma::mat densityMatrix;
    arma::mat C;

    Basis *basis = nullptr;

    void setCMatrix();
    void updateDensityMatrix();
public:
    HartreeFock();
    HartreeFock(int NElectrons, int singleParticleStates);
    ~HartreeFock();
    void initializeHF(int NElectrons, int singleParticleStates, Basis *newBasis);
    int runHF(int maxHFIterations);
    void updateHFMatrix(arma::mat &HFMatrix);
    double calculateInnerHFMatrixElement(int alpha, int alpha_ml, int alphaSpin, int beta, int beta_ml, int betaSpin);
    void writeToFile(arma::vec eigVals, arma::mat eigVecs);
    void getHFEnergy();

    // Printers
    void printHFMatrix(arma::mat HFMatrix);
    void printSPEnergies(arma::vec singleParticleEnergies);

    // Setters
    void setConvergence(double newLambda) { lambda = newLambda; }
    void setInteractionMatrix(double * newInteractionMatrix);
};

#endif // HARTREEFOCK_H
