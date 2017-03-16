#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <cmath>
#include <armadillo>
#include "basis.h"
#include "hermitepolynomials.h"
#include "Coulomb_Functions.h"

class HartreeFock
{
private:
    int N_SPS;
    int N_Electrons;
    int HFCounter = 0;
    double lambda = 1e-5;
    double sqrtOmega = 1;
    double * interactionMatrix;
    arma::vec SPS_Energies;
    arma::mat densityMatrix;
    arma::mat C;

    Basis *basis = nullptr;

    void initializeCMatrix();
    void updateDensityMatrix();
public:
    HartreeFock();
    HartreeFock(int NElectrons, int singleParticleStates);
    ~HartreeFock();
    void initializeHF(int NElectrons, int singleParticleStates, Basis *newBasis);
    int runHF(int maxHFIterations);
    void updateHFMatrix(arma::mat &HFMatrix);
    double calculateInnerHFMatrixElement(int alpha, int alpha_ml, int alphaSpin, int beta, int beta_ml, int betaSpin);
    void writeToFile();
    void getHFEnergy(double &HFEnergyResults, int &HFIterationsResults);

    // Printers
    void printHFMatrix(arma::mat HFMatrix);
    void printSPEnergies();

    // Setters
    void setConvergence(double newLambda) { lambda = newLambda; }
    void setInteractionMatrix(double * newInteractionMatrix);
    void setOmega(double newOmega) { sqrtOmega = sqrt(newOmega); std::cout<<sqrtOmega<<std::endl; }
};

#endif // HARTREEFOCK_H
