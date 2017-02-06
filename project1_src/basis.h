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

    void initializeBasis(int ECutoff);

    int getShellDegeneracy(int shellNumber);
    int getTotalDegeneracy();
    int getTotalParticleNumber();
    int getShellParticles(int shellNumber);
    double getTotalEnergy();
    double getShellEnergy(int shellNumber);

    SingleState * getState(int state);

    // getDegeneracyOfShell
    // getTotalEnergy
    // getParticleNumber
    // getShellParticles
};

#endif // BASIS_H
