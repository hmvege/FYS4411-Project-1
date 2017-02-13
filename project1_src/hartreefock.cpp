#include "hartreefock.h"

HartreeFock::HartreeFock(int newNElectrons, int newMaxShell)
{    
    N_Electrons = newNElectrons;
    maxShell = newMaxShell;
}

HartreeFock::~HartreeFock()
{

}

void HartreeFock::setInteractionMatrix(double * newInteractionMatrix)
{
    interactionMatrix = newInteractionMatrix;
}
