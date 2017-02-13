#ifndef UNITTESTS_H
#define UNITTESTS_H

int runHermiteTimer(int N, double X);
int checkCorrectHermitePolynoms(int maxHermite);
void testIntegratorClass(int NIntPoints, double eps=1e-10);
void testIntegratorFunction(int NIntPoints, double eps=1e-10);

#endif // UNITTESTS_H
