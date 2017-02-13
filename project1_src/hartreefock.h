#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include "basis.h"
#include "hermitepolynomials.h"

class HartreeFock
{
private:
    int N_Electrons;
    int maxShell;
    double * interactionMatrix;

    HermitePolynomials hermite;
    Basis basis;
public:
    HartreeFock(int N_Electrons, int newMaxShell);
    ~HartreeFock();

    void setInteractionMatrix(double * newInteractionMatrix);


};

#endif // HARTREEFOCK_H
