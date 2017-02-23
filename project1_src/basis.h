#ifndef BASIS_H
#define BASIS_H
#include <vector>
#include <singlestate.h>

using std::vector;

class Basis
{
private:
    vector<SingleState *> statesArray;
public:
    Basis();
    ~Basis();

    double omega;

    void initializeBasis(int ECutoff, double initOmega);

    int getShellDegeneracy(int shellNumber);
    int getTotalDegeneracy();
    int getTotalParticleNumber();
    int getShellParticles(int shellNumber);
    double getOmega();
    double getTotalEnergy();
    double getShellEnergy(int shellNumber);
    void printBasisSize();
    void printBasis(bool polar=false);

    SingleState * getState(int state);

    // getDegeneracyOfShell
    // getTotalEnergy
    // getParticleNumber
    // getShellParticles
};

#endif // BASIS_H
