#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"rkf45.h"

/*
Example of use.
*/

void testFunc(double t, double *x, double *dx, PARAMETERS* p)
{

    dx[0] = p->doubleVector[0]*sqrt(t);
    dx[1] = p->doubleVector[1]*x[0];
    dx[2] = p->doubleVector[2]*x[1];    

}

int main()
{

    double x0[3] = {0.0, 0.0, 0.0};
    FILE* ff = fopen("output.csv", "w");
    PHASE* p = phaseInitialize();
    p->parameters->doubleVector = (double*)malloc(3*sizeof(double));
    
    p->tSim = 1.;
    p->nStates = 3;    
    p->parameters->doubleVector[0] = 3./2;
    p->parameters->doubleVector[1] = 5./2;
    p->parameters->doubleVector[2] = 7./2;        

    phaseIntegrate(p, x0, testFunc, ff);
    
    fclose(ff);

    return 0;

}

