#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <armadillo>
#include "basis.h"
#include "hermitepolynomials.h"

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
    void writeToFile(arma::vec eigVals, arma::mat eigVecs);
    void getEnergies();

    // Setters
    void setConvergence(double newLambda) { lambda = newLambda; }
    void setInteractionMatrix(double * newInteractionMatrix);
};

#endif // HARTREEFOCK_H
