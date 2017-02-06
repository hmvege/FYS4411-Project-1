#include "basis.h"
#include <vector>

using std::vector;
using std::cout;
using std::endl;

Basis::Basis()
{

}

Basis::~Basis()
{
    // Destructor
    for (SingleState * state : statesArray)
    {
        delete state;
    }
    statesArray.clear();
}

void Basis::initializeBasis(int ECutoff)
{
    // Retrieving all possible states up until a cutoff-shell

    for (int i = 0; i < ECutoff; i++) // Running over variations of nx
    {
        for (int j = 0; j < ECutoff; j++) // Running over variations of ny
        {
            for (double k = -0.5; k < 1; k++) // Running over the two possible spin configurations
            {
                if (i+j < ECutoff) // If-test to prevent counting outside of the shell
                {
                SingleState * state = new SingleState();
                state->set(i, j, k);
                statesArray.push_back(state);
                }
            }
        }
    }
}

SingleState *Basis::getState(int state)
{
//    cout << statesArray << endl;
    return statesArray[state];
}

int Basis::getShellDegeneracy(int shellNumber)
{
    int degeneracy = 0;
    for (SingleState * state : statesArray)
    {
        if ((state->getN_x() + state->getN_y()) == (shellNumber - 1))
        {
            degeneracy++;
        }
    }
    return degeneracy;
}

int Basis::getTotalDegeneracy()
{
    return statesArray.size();
}


double Basis::getShellEnergy(int shellNumber)
{
    int energy = 0;
    for (SingleState * state : statesArray)
    {
        if ((state->getN_x() + state->getN_y()) == (shellNumber - 1))
        {
            energy += state->getEnergy();
        }
    }
    return energy;
}

double Basis::getTotalEnergy()
{
    double energy = 0;

    for (SingleState * state : statesArray)
    {
        energy += state->getEnergy();
    }

    return energy;
}

int Basis::getTotalParticleNumber()
{
    return statesArray.size();
}

int Basis::getShellParticles(int shellNumber)
{
    int particleNumber = 0;

    for (SingleState * state : statesArray)
    {
        if ((state->getN_x() + state->getN_y()) == (shellNumber - 1))
        {
            particleNumber++;
        }
    }
    return particleNumber;
}
