#include "singlestate.h"

using std::cout;
using std::endl;

SingleState::SingleState()
{
    int n_x = 0;
    int n_y = 0;
    double spin = 0;
    double energy = 0;
}


SingleState::SingleState(int initN_x, int initN_y, double initSpin)
{
    int n_x = initN_x;
    int n_y = initN_y;
    double spin = initSpin;
}

void SingleState::printSystem()
{
    cout << "n_x = " << n_x << ", n_y = " << n_y << ", spin = " << spin << ", energy = " << getEnergy() << endl;
}

void SingleState::set(int newN_x, int newN_y, double newSpin)
{
    setN_x(newN_x);
    setN_y(newN_y);
    setSpin(newSpin);
}
