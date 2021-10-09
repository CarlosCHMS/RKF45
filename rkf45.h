#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct{

    /*
    This struct permits pass parametric information to
    the function thats going to be integrated.
    
    It is possible pass double, int and char variables.
    
    */

    double* doubleVector;
    int* intVector;
    char* string;    
    
} PARAMETERS;

typedef struct{

    // Number of states
    int nStates;
    
    // Duration of the simulation
    double tSim;
    
    // Time step
    double dt;
    
    // Maximal time step
    double dtMax;    
    
    //Flag to save each step in a .csv file
    int save;
    
    // Final states
    double* xFinal;
    
    // Final time
    double tFinal;
    
    //Local error tolerance
    double tol;
    
    // Funtion paramters
    PARAMETERS* parameters;
    
} PHASE;

PHASE* phaseInitialize();

void phasePropagateRKF45(PHASE *phase, double* t, double *x, void (*func)(double, double*, double*, PARAMETERS*));

void phaseIntegrate(PHASE *phase, double *x0, void (*func)(double, double*, double*, PARAMETERS*) , FILE* ff);

