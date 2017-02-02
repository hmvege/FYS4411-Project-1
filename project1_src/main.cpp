#include <iostream>
#include <vector>
#include "singlestate.h"
#include "basis.h"
#include "hermitepolynomials.h"

#include "unittests.h"
//#include <iomanip>
//#include <time.h>

using namespace std;

void *getStatesArray();

int main(int nargs, char *args[])
{
    Basis basis;
    basis.initializeBasis(10);

    cout <<"Shell 4:" << endl;
    cout <<"ShellPar: "<< basis.getShellParticles(4) << endl;
    cout <<"ShellDeg: "<< basis.getShellDegeneracy(4) << endl;
    cout <<"ShellEne: "<< basis.getShellEnergy(4) << endl;

    cout <<"Total:" << endl;
    cout <<"TotalPar: "<< basis.getTotalParticleNumber() << endl;
    cout <<"TotalDeg: "<< basis.getTotalDegeneracy() << endl;
    cout <<"TotalEne: "<< basis.getTotalEnergy() << endl;


    // Tests
    runHermiteTimer(1e3, 5.5);
//    checkCorrectHermitePolynoms(6);

    return 0;
}
