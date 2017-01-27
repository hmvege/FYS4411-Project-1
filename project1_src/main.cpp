#include "singlestate.h"
#include <iostream>
#include <vector>

using namespace std;

void *getStatesArray();

int main(int nargs, char *args[])
{
//    SingleState state1;
//    state1.setN_x(0);
//    state1.setN_y(0);
//    state1.setSpin(-0.5);
//    state1.setEnergy(1);
//    state1.printSystem();

    getStatesArray();

    return 0;
}

void * getStatesArray()
{
    vector<SingleState *> statesArray;
    int maxShell = 10;

    for (int i = 0; i < maxShell; i++) // Running over variations of nx
    {
        for (int j = 0; j < maxShell; j++) // Running over variations of ny
        {
            for (double k = -0.5; k < 1; k++) // Running over the two possible spin configurations
            {
                if (i+j < maxShell)
                {
                SingleState * state = new SingleState();
                state->set(i, j, k);
                statesArray.push_back(state);
                }
            }
        }
    }


    cout << statesArray.size() << endl;
//    for (SingleState * state : statesArray)
//    {
//        state->printSystem();
//    }

    return 0;
}
