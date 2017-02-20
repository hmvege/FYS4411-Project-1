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

    // Setters
    void setInteractionMatrix(double * newInteractionMatrix);
//    void setN_SPS(int N) { N_SPS = N; }
//    void setN_Electrons(int N) { N_Electrons = N; }
};

#endif // HARTREEFOCK_H
