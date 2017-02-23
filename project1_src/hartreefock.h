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

    double lambda = 1e-5; // TODO: Add this as a user setting in QuantumDot class

    double * interactionMatrix;
//    double ** densityMatrix;
//    double ** C;

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

    // Setters
    void setConvergence(double newLambda) { lambda = newLambda; }
    void setInteractionMatrix(double * newInteractionMatrix);
//    void setN_SPS(int N) { N_SPS = N; }
//    void setN_Electrons(int N) { N_Electrons = N; }
};

#endif // HARTREEFOCK_H
