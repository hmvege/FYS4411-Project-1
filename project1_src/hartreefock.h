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
    int maxHFIterations = 100;

    double * interactionMatrix;
//    double ** densityMatrix;
//    double ** C;

    arma::mat densityMatrix;
    arma::mat C;


    HermitePolynomials hermite;
    Basis basis;
public:
    HartreeFock();
    HartreeFock(int NElectrons, int singleParticleStates);
    ~HartreeFock();
    void initializeHF(int NElectrons, int singleParticleStates);
    int runHF();

    // Setters
    void setInteractionMatrix(double * newInteractionMatrix);
    void setMaxHFIterations(int N) { maxHFIterations = N; }
//    void setN_SPS(int N) { N_SPS = N; }
//    void setN_Electrons(int N) { N_Electrons = N; }
};

#endif // HARTREEFOCK_H
