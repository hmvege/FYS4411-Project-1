#ifndef SINGLESTATE_H
#define SINGLESTATE_H

#include <iostream>

class SingleState
{
public:
    int n_x;
    int n_y;
    double spin;

    SingleState();
    SingleState(int initN_x, int initN_y, double initSpin);

    // Setters
    void setN_x(int newN_x) { n_x = newN_x; }
    void setN_y(int newN_y) { n_y = newN_y; }
    void setSpin(double newSpin) { spin = newSpin; }
    void set(int newN_x, int newN_y, double newSpin);

    // Getters
    double getN_x() { return n_x; }
    double getN_y() { return n_y; }
    double getSpin() { return spin; }
    double getEnergy() { return (n_x + n_y + 1); }

    // Printers
    void printSystem();
};

#endif // SINGLESTATE_H
