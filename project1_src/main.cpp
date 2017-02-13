#include <iostream>
#include <cmath>
#include <vector>
#include "singlestate.h"
#include "basis.h"
#include "hermitepolynomials.h"
#include "functions.h"

#include "quantumdot.h"

#include "hartreefock.h"

#include <armadillo>
#include "unittests.h"
#include <iomanip>
//#include <time.h>

using namespace std;

double potentialV(double x1, double x2, double y1, double y2);
void printMatrix(double ** A, int N);
void setMatrixZero(double ** A, int N);

int main(int nargs, char *args[])
{
    int N = 2; // Number of electrons, should be magic number: 2, 6, 12, 20
    int maxShell = 4;

//    testIntegratorClass(2); // Works
//    testIntegratorFunction(2); // Works

    HermitePolynomials hermite;

    Basis basis;
    basis.initializeBasis(maxShell);
    int N_SPS = basis.getTotalParticleNumber();

    // HF setup
    double * interactionMatrix = new double[(int) pow(N_SPS,4)];
    double ** densityMatrix     = new double * [N_SPS];
    double ** C                 = new double * [N_SPS];

    for (int i = 0; i < N_SPS; i++)
    {
        C[i]                    = new double[N_SPS];
        densityMatrix[i]        = new double[N_SPS];
    }

    // Initializing C as a diagonal matrix(other initializations exists)
    for (int i = 0; i < N_SPS; i++)
    {
        C[i][i] = 1;
    }

    // Setting up the density matrix
    for (int gamma = 0; gamma < N_SPS; gamma++)
    {
        for (int delta = 0; delta < N_SPS; delta++)
        {
            double sum = 0;
            for (int i = 0; i < N; i++)
            {
                sum += C[i][gamma]*C[i][delta];
            }
            densityMatrix[gamma][delta] = sum;
        }
    }

    // Setting up interaction matrix
    for (int alpha = 0; alpha < N_SPS; alpha++)
    {
        double (*alphaHerm_x)(double) = hermite.getPolynom(basis.getState(alpha)->getN_x());
        double (*alphaHerm_y)(double) = hermite.getPolynom(basis.getState(alpha)->getN_y());
        double alphaConst = basis.getState(alpha)->normConst();
        for (int beta = 0; beta < N_SPS; beta++)
        {
            double (*betaHerm_x)(double) = hermite.getPolynom(basis.getState(beta)->getN_x());
            double (*betaHerm_y)(double) = hermite.getPolynom(basis.getState(beta)->getN_y());
            double betaConst = basis.getState(beta)->normConst();
            for (int gamma = 0; gamma < N_SPS; gamma++)
            {
                double (*gammaHerm_x)(double) = hermite.getPolynom(basis.getState(gamma)->getN_x());
                double (*gammaHerm_y)(double) = hermite.getPolynom(basis.getState(gamma)->getN_y());
                double gammaConst = basis.getState(gamma)->normConst();
                for (int delta = 0; delta < N_SPS; delta++)
                {
                    double (*deltaHerm_x)(double) = hermite.getPolynom(basis.getState(delta)->getN_x());
                    double (*deltaHerm_y)(double) = hermite.getPolynom(basis.getState(delta)->getN_y());
                    double deltaConst = basis.getState(delta)->normConst();

                    interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)] = solveGaussianHermiteQuadrature(maxShell,
                                                                                                                alphaHerm_x,
                                                                                                                alphaHerm_y,
                                                                                                                betaHerm_x,
                                                                                                                betaHerm_y,
                                                                                                                gammaHerm_x,
                                                                                                                gammaHerm_y,
                                                                                                                deltaHerm_x,
                                                                                                                deltaHerm_y,
                                                                                                                potentialV);
                    interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)] *= alphaConst*betaConst*gammaConst*deltaConst;
                    // state_alpha*alpha_norm * etc..
                }
            }
        }
    }

//    cout << interactionMatrix[0] << endl;

    // For loop will run till max HF iteration is reached, or we get a convergence for epsilon(NOT IMPLEMENTED YET)
    int maxHFiterations = 100;
    for (int HFiteration = 0; HFiteration < maxHFiterations; HFiteration++)
    {
        // Setting up HFMatrix
        double ** HFMatrix = new double * [N_SPS];
        for (int i = 0; i < N_SPS; i++) {
            HFMatrix[i] = new double[N_SPS];
        }
        for (int i = 0; i < N_SPS; i++) {
            for (int j = 0; j < N_SPS; j++) {
                HFMatrix[i][j] = 0;
            }
        }

        // retrieving HF elements
        for (int alpha = 0; alpha < N_SPS; alpha++)
        {
            for (int beta = 0; beta < N_SPS; beta++)
            {
                double HFElement = 0;
                for (int gamma = 0; gamma < N_SPS; gamma++)
                {
                    for (int delta = 0; delta < N_SPS; delta++)
                    {
                        HFElement += C[gamma][delta] * interactionMatrix[index(alpha, beta, gamma, delta, N_SPS)];
                    }
                }
                HFMatrix[alpha][beta] = HFElement;
                if (alpha == beta) // Delta function
                {
                    HFMatrix[alpha][beta] += basis.getState(alpha)->getEnergy();
                }
            }
        }

        // Finding eigenvalues
        arma::mat A = arma::zeros<arma::mat>(N_SPS,N_SPS);
        // When if we have convergence, rather than checking everyone, find the max element

        // De-allocating memory
        for (int i = 0; i < N_SPS; i++) {
            delete [] HFMatrix[i];
        }
        delete [] HFMatrix;
    }





//    // TEST: Setting up basis
//    Basis basis;
//    basis.initializeBasis(10);

//    // TEST: Checking if getters are working
//    cout <<"Shell 4:" << endl;
//    cout <<"ShellPar: "<< basis.getShellParticles(4) << endl;
//    cout <<"ShellDeg: "<< basis.getShellDegeneracy(4) << endl;
//    cout <<"ShellEne: "<< basis.getShellEnergy(4) << endl;

//    cout <<"Total:" << endl;
//    cout <<"TotalPar: "<< basis.getTotalParticleNumber() << endl;
//    cout <<"TotalDeg: "<< basis.getTotalDegeneracy() << endl;
//    cout <<"TotalEne: "<< basis.getTotalEnergy() << endl;

//    // TEST: Comparing time usage and number of correct hermite polynomials
//    runHermiteTimer(1e6, 5.5);
//    checkCorrectHermitePolynoms(6);


    // Deallocation
    for (int i = 0; i < N; i++)
    {
        delete [] densityMatrix[i];
//        delete [] C[i];
    }
    delete [] interactionMatrix;
    delete [] densityMatrix;
    delete [] C;

    return 0;
}

void printMatrix(double ** A, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << A[i][j] << " ";
//            cout << A[i*N+j] << " ";
        }
        cout << endl;
    }
}

double potentialV(double x1, double x2, double y1, double y2)
{
    return 1.0/sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

void setMatrixZero(double ** A, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = 0;
        }
    }
}
