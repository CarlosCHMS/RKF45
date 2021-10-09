#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"rkf45.h"


/*
Runge–Kutta–Fehlberg method obtained from:
SOME EXPERIMENT&L RESULTS CONCERNING
THE ERROR PROPAGATION IN RUNGE-KUTTA
TYPE INTEGRATION FORMULAS
*/


PHASE* phaseInitialize()
{

    PHASE* phase = malloc(sizeof(PHASE));  

    phase->xFinal = malloc(phase->nStates*sizeof(double));
    
    phase->parameters = (PARAMETERS*)malloc(sizeof(PARAMETERS));
    
    phase->save = 1;
    
    phase->tol = 1e-8;
    
    return phase;

}


void phasePropagateRKF45(PHASE *phase, double* t, double *x, void (*func)(double, double*, double*, PARAMETERS*))
{

    int jj;
    int cont;
    
    double tAux;
    double TE;
    double h;
    double TEaux;
    
    double k1[phase->nStates];
    double k2[phase->nStates];
    double k3[phase->nStates];
    double k4[phase->nStates];
    double k5[phase->nStates];
    double k6[phase->nStates];        
    double xAux[phase->nStates];

    h = phase->dt;

    cont = 1;
    while(cont)
    {

        //Step 1    
        func(*t, x, k1, phase->parameters);

        //Step 2
        tAux = *t + (1./4)*h;
        for(jj=0; jj<phase->nStates; jj++)
        {       
            xAux[jj] = x[jj] + (1./4)*h*k1[jj];
                    
        }        
        func(tAux, xAux, k2, phase->parameters);    
        
        //step3
        tAux = *t + (3./8)*h;
        for(jj=0; jj<phase->nStates; jj++)
        {       
            xAux[jj] = x[jj] + (3./32)*h*k1[jj] + (9./32)*h*k2[jj];
            
        }    
        func(tAux, xAux, k3, phase->parameters);        
        
        //step4
        tAux = *t + (12./13)*h;
        for(jj=0; jj<phase->nStates; jj++)
        {       
            xAux[jj] = x[jj] + (1932./2197)*h*k1[jj] + (-7200./2197)*h*k2[jj] + (7296./2197)*h*k3[jj];
            
        }    
        func(tAux, xAux, k4, phase->parameters);    

        //step5
        tAux = *t + h;
        for(jj=0; jj<phase->nStates; jj++)
        {       
            xAux[jj] = x[jj] + (439./216)*h*k1[jj] + (-8.)*h*k2[jj] + (3680./513)*h*k3[jj] + (-845./4104)*h*k4[jj];
            
        }    
        func(tAux, xAux, k5, phase->parameters);    

        //step6
        tAux = *t + (1./2)*h;
        for(jj=0; jj<phase->nStates; jj++)
        {       
            xAux[jj] = x[jj] + (-8./27)*h*k1[jj] + (2.)*h*k2[jj] + (-3544./2565)*h*k3[jj] + (1859./4104)*h*k4[jj] + (-11./40)*h*k5[jj];
            
        }    
        func(tAux, xAux, k6, phase->parameters);    
       
        //Error
        TE = 0.;
        for(jj=0; jj<phase->nStates; jj++)
        {       
            TEaux = h*((-1./360)*k1[jj] + (128./4275)*k3[jj] + (2197./75240)*k4[jj] + (-1./50)*k5[jj] + (-2./55)*k6[jj]);
            if(TE < fabs(TEaux))
            {
                TE = fabs(TEaux);
            }
        }
        
        if(TE > phase->tol)
        {
            h = 0.9*h*pow(phase->tol/fabs(TE), 0.2);
        }
        else
        {
            cont = 0;
        }
    
    }
        
    //Propagation
    *t += h;
    for(jj=0; jj<phase->nStates; jj++)
    {       
        x[jj] += h*((16./135)*k1[jj] + (6656./12825)*k3[jj] + (28561./56430)*k4[jj] + (-9./50)*k5[jj] + (2./55)*k6[jj]);
        
    }
    
    //Next step
    if(TE > 0)
    {
        phase->dt = 0.9*h*pow(phase->tol/fabs(TE), 0.2);
    }
    
    //Max step size
    if(phase->dt > phase->dtMax)
    {
        phase->dt = phase->dtMax;    
    }
    
}


void phaseIntegrate(PHASE *phase, double *x0, void (*func)(double, double*, double*, PARAMETERS*) , FILE* ff)
{

    int ii, jj;
    
    double t = 0.0;
    double x[phase->nStates];
    
    phase->dtMax = phase->tSim/10;
    phase->dt = phase->dtMax;

    for(jj=0; jj<phase->nStates; jj++)
    {
        x[jj] = x0[jj];    
    }

    //Save initial values
    if(phase->save)
    {  
        fprintf(ff, "%+e,", t);        
                    
        for(jj=0; jj<phase->nStates; jj++)
        {
            fprintf(ff, "%+e,", x[jj]);             
        }            
        fprintf(ff, "\n");
    }

    while(t < phase->tSim)
    {        
        //Final step
        if(t + phase->dt > phase->tSim)
        {
            phase->dt = phase->tSim - t;            
        }
        
        //Propagator
        phasePropagateRKF45(phase, &t, x, func);
        
        //Save
        if(phase->save)
        {           
            fprintf(ff, "%+e,", t);        
                        
            for(jj=0; jj<phase->nStates; jj++)
            {
                fprintf(ff, "%+e,", x[jj]);             
            }            
            fprintf(ff, "\n");
        }
    }
    
    //Final state and time
    for(jj=0; jj<phase->nStates; jj++)
    {
        phase->xFinal[jj] = x[jj];    
    }    

    phase->tFinal = t;

}

