#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H

#include "hermitepolynomials.h"
#include "basis.h"

class quantumDot
{
private:
    int N_Electrons;
    int maxShell;
    int N_SPS;
    double * interactionMatrix;
    double (*V)(double x1, double x2, double y1, double y2);
    HermitePolynomials hermite;
    Basis basis;
public:
    quantumDot(int newNElectrons, int newMaxShell);

    void setPotential(double (*newV)(double x1, double x2, double y1, double y2));
    void setupInteractionMatrix();
    void runHartreeFock();

};

#endif // QUANTUMDOT_H
