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
public:
    quantumDot(int newNElectrons, int newMaxShell, double newOmega);
    ~quantumDot() { }

    void setHFLambda(double newLambda) { HF.setConvergence(newLambda); }
    void setPotential(double (*newV)(double x1, double x2, double y1, double y2));
    void setupInteractionMatrix(int integrationPoints);
    void setupInteractionMatrixPolar();
    void setupInteractionMatrixFromFile(const std::string& filename);
    void runHartreeFock(int maxHFIteration);
    void printInteractionMatrix(int NPrintPoints);
};

#endif // QUANTUMDOT_H
