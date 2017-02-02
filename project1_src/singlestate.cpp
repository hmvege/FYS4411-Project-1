#include "singlestate.h"

using std::cout;
using std::endl;

SingleState::SingleState()
{
    n_x = 0;
    n_y = 0;
    spin = 0;
}


SingleState::SingleState(int initN_x, int initN_y, double initSpin)
{
    n_x = initN_x;
    n_y = initN_y;
    spin = initSpin;
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
