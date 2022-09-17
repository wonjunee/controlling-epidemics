#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <iostream>


void intialize_rho0(double* rho, const int n1, const int n2, const int nt){

	double sum=0;
	double vmax=0;
    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;

            /* Exp1 */
            rho[i*n1+j] =  exp(-5*pow(x-0.5,2)-5*pow(y-0.5,2))*0.5;
        }
    }


    for(int n=1;n<nt;++n){
        for(int i=0;i<n1*n2;++i){
            rho[n*n1*n2+i]=rho[i];
        }
    }
}

void intialize_rho1(double* rho, const int n1, const int n2, const int nt){

    double sum=0;
    double vmax = 0;

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){

            double x=1.0*(j+0.5)/n1;
            double y=1.0*(i+0.5)/n2;

            // experiment 1
            rho[i*n1+j] =  0.2* exp(-30*pow(x-0.64,2)-30*pow(y-0.64,2));
            rho[i*n1+j] += 0.2* exp(-30*pow(x-0.45,2)-30*pow(y-0.45,2));
        }
    }

    for(int n=1;n<nt;++n){
        for(int i=0;i<n1*n2;++i){
            rho[n*n1*n2+i]=rho[i];
        }
    }
}

void intialize_rho2(double* rho, const int n1, const int n2, const int nt){

    for(int i=0;i<n2;++i){
        for(int j=0;j<n1;++j){
            rho[i*n1+j] = 1e-2;
        }
    }

    for(int n=1;n<nt;++n){
        for(int i=0;i<n1*n2;++i){
            rho[n*n1*n2+i]=rho[i];
        }
    }
}

#endif