#ifndef SINGLESTATE_H
#define SINGLESTATE_H

#include <iostream>
#include <cmath>
#include <math.h>
#include "hermitepolynomials.h"
#include "functions.h"

class SingleState
{
private:
    int n_x;
    int n_y;
    int n;
    int m;
    double spin;
    double omega;

    void setPolar(int newN_x, int newN_y);
public:
    SingleState();
    SingleState(int initN_x, int initN_y, double initSpin, double initEnergy);

    // Setters for Cartesian quantum numbers
    void setN_x(int newN_x) { n_x = newN_x; }
    void setN_y(int newN_y) { n_y = newN_y; }
    void setSpin(double newSpin) { spin = newSpin; }
    void setOmega(double newOmega) { omega = newOmega; }
    void set(int newN_x, int newN_y, double newSpin, double omega);

    // Getters
    double getN_x() { return n_x; }
    double getN_y() { return n_y; }
    double getN() { return n; }
    double getM() { return m; }
    double getSpin() { return spin; }
    double getEnergy() { return omega*(n_x + n_y + 1); }

//    double normConst() { return 1./sqrt(pow(2,((double) n_x) + ((double) n_y))) * ((double) factorial(n_x)) * ((double) factorial(n_y)); }
    double normConst() { return 1./sqrt(pow(2,((double) n_x) + ((double) n_y)) * ((double) factorial(n_x)) * ((double) factorial(n_y)) * M_PI); }

    // Printers
    void printSystem(bool polar);
};

#endif // SINGLESTATE_H
