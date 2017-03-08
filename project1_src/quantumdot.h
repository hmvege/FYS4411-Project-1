#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H

#include "hermitepolynomials.h"
#include "basis.h"
#include "hartreefock.h"

class quantumDot
{
private:
    int N_Electrons;
    int maxShell;
    int N_SPS;
    int interactionMatrixLength;
    double omega;
    double * interactionMatrix;
    double (*V)(double x1, double x2, double y1, double y2);

    HermitePolynomials hermite;
    Basis basis;
    HartreeFock HF;

    // Variables for storing Hartree-Fock results
    double HFEnergyResults = 0;
    int HFIterationsResults = 0;
public:
    quantumDot(int newNElectrons, int newMaxShell, double newOmega);
    ~quantumDot() { }

    double runHartreeFock(int maxHFIteration);
    void antiSymmetrizeMatrix();
    void setupInteractionMatrixFromFile(const std::string& filename);
    void setupInteractionMatrix(int integrationPoints);
    void setupInteractionMatrixPolar();
    void printInteractionMatrix(int NPrintPoints);
    void storeResults(const std::string& filename);

    // Setters
    void setHFLambda(double newLambda) { HF.setConvergence(newLambda); }
    void setPotential(double (*newV)(double x1, double x2, double y1, double y2));

    // Getters
    int getInteractionMatrixLength() { return interactionMatrixLength; }
    Basis getBasis() { return basis; }
};

#endif // QUANTUMDOT_H
