#ifndef METHOD_H
#define METHOD_H
#include <iostream>

#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include <omp.h> // parallelization
#include <fftw3.h>
#include "poisson_solver_3d.h"
#include "helper.h"
#include <memory>

class Method{
public:
    int n1;
    int n2;
    int nt;

    // For SIR model
    double beta;
    double gamma;

    int im;
    int ip;
    int jm;
    int jp;

    double lambda_rho; // parameter for strong convexity of rho
    double lambda_m;   // parameter for strong convexity of m

    double tau[3];     // stepsize for primal variables
    double sigma[3];   // stepsize for dual variables

    int max_iteration; // maximum iteration of the PDHG
    double tolerance;  // tolerance for the error (stopping condition)

    double* mx[3];
    double* my[3];
    double* phi[3];
    double* phitmps[3];
    double* rhotmps[3];

    double* rhont0tmps[3];

    double* psi[3];

    double etalist[3];
    double alphalist[3];

    unique_ptr<poisson_solver>     fftps[3];
    unique_ptr<poisson_solver_DST> fftpsDST;

    Method(){
        for(int i=0;i<3;++i){
            phitmps[i]=nullptr;
            rhotmps[i]=nullptr;
            rhont0tmps[i]=nullptr;
            mx[i]=nullptr;
            my[i]=nullptr;
            phi[i]=nullptr;
            psi[i]=nullptr;
        }
    }

    Method(int n1, int n2, int nt, double tau, double sigma, int max_iteration, double tolerance, 
        double alphalist[3], double* etalist, double beta, double gamma)
    :   n1(n1), n2(n2), nt(nt), max_iteration(max_iteration), tolerance(tolerance),
        beta(beta), gamma(gamma)

    {

        lambda_rho = 0.1;  // lambda_rho strong convexity
        lambda_m   = 1000; // lambda_rho strong convexity

        this->tau[0] = tau;
        this->tau[1] = tau;
        this->tau[2] = tau;

        this->sigma[0] = sigma;
        this->sigma[1] = sigma;
        this->sigma[2] = sigma;

        this->alphalist[0] = alphalist[0];
        this->alphalist[1] = alphalist[1];
        this->alphalist[2] = alphalist[2];

        for(int i=0;i<3;++i){
            mx[i]         = new double[n1*n2*nt];
            my[i]         = new double[n1*n2*nt];
            phi[i]        = new double[n1*n2*nt];
            phitmps[i]    = new double[n1*n2*nt];
            rhotmps[i]    = new double[n1*n2*nt];
            rhont0tmps[i] = new double[n1*n2];
            psi[i]        = new double[n1*n2];
        }

        clock_t t;
        t = clock();
            
        double eta = etalist[0];

        // initialize poisson solvers
        for(size_t k=0;k<3;++k){
            fftps[k] = unique_ptr<poisson_solver>(new poisson_solver(n1,n2,nt,eta));
        }
        fftpsDST = unique_ptr<poisson_solver_DST>(new poisson_solver_DST(n1,n2,nt,eta));

        t = clock() - t;
        printf ("\nCPU time for setting up FFT: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
    }


    ~Method(){
        for(int k=0;k<3;++k){
            delete[] mx[k];
            delete[] my[k];
            delete[] rhotmps[k];
            delete[] rhont0tmps[k];
            delete[] phitmps[k];
            delete[] phi[k];
            delete[] psi[k];
        }
    }

    void setup_indices(int& im, int& ip, int& jm, int& jp, const int i, const int j){
        im = fmax(0,i-1);
        ip = fmin(n2-1,i+1);

        jm = fmax(0,j-1);
        jp = fmin(n1-1,j+1);
    }

    void perform_upwind_scheme(double& muxp, double& muxm, double& muyp, double& muym, const double* phi, const int i, const int j){

        setup_indices(im,ip,jm,jp,i,j);

        muxp = 1.0 * n1 * (phi[i*n1+jp] - phi[i*n1+j]);
        muxm = 1.0 * n1 * (phi[i*n1+j] - phi[i*n1+jm]);
        muyp = 1.0 * n2 * (phi[ip*n1+j] - phi[i*n1+j]);
        muym = 1.0 * n2 * (phi[i*n1+j] - phi[im*n1+j]);
    }

    // centered difference
    void update_m(double* mx, double* my, const double* rho, const double* phi, const int SIR){

        int im,ip,jm,jp;

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;

                    setup_indices(im,ip,jm,jp,i,j);

                    double nablaxphi = 1.0*n1*(phi[n*n1*n2+i*n1+jp]-phi[ind]);
                    double nablayphi = 1.0*n2*(phi[n*n1*n2+ip*n1+j]-phi[ind]);
                    
                    double rhovalx=rho[ind];
                    double rhovaly=rho[ind];

                    // mx[n*n1*n2+i*n1+j] = rhovalx/(tau[SIR] * alphalist[SIR] + rhovalx*(1+tau[SIR]+lambda_rho)) * (mx[n*n1*n2+i*n1+j] - tau[SIR] * nablaxphi);
                    // my[n*n1*n2+i*n1+j] = rhovaly/(tau[SIR] * alphalist[SIR] + rhovaly*(1+tau[SIR]+lambda_rho)) * (my[n*n1*n2+i*n1+j] - tau[SIR] * nablayphi);
                    mx[ind] = (1.0/tau[SIR] * mx[ind] - nablaxphi)/(1.0/rhovalx + lambda_m + 1.0/tau[SIR]);
                    my[ind] = (1.0/tau[SIR] * my[ind] - nablayphi)/(1.0/rhovaly + lambda_m + 1.0/tau[SIR]);
                }
            }   
        }   
    }

    /* 
        Update rho
    */

    void calculate_rho_related(double& mvalue, double& Dtphi, double& Deltaphi, const int n, const int i, const int j, const double* mx, const double* my, const double* phi){

        double mxvalue, myvalue;
        int im,ip,jm,jp;
        setup_indices(im,ip,jm,jp,i,j);

        int ind = n*n1*n2+i*n1+j;

        mxvalue = mx[ind];
        myvalue = my[ind];
        
        mvalue = sqrt(mxvalue*mxvalue + myvalue*myvalue);

        if(n<nt-1) Dtphi=1.0*nt*(phi[(n+1)*n1*n2+i*n1+j]-phi[ind]);
        else       Dtphi=0;
        Deltaphi = - n1*n1 * (-phi[n*n1*n2+i*n1+jm]+2*phi[ind]-phi[n*n1*n2+i*n1+jp])
                   - n2*n2 * (-phi[n*n1*n2+im*n1+j]+2*phi[ind]-phi[n*n1*n2+ip*n1+j]);

    }

    void update_rho0(double* rho0,const double* rho1, const double* rho2, const double* rho3,const double* mx,const double* my,const double* obstacle){

    	double newrhovalue = -1;
		double tauval = tau[0];
		double betaval  = beta;

        // memcpy(&rho0[0], rhont0tmps[0], n1*n2*sizeof(double));

		for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;
                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;
                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[0]);
                    double aval = 0;
                    double cval = -0.5*alphalist[0]*mvalue*mvalue / (lambda_rho + 1.0/tauval);

                    if(n==nt-1){
                        aval = 1.0/(lambda_rho + 1.0/tauval) * 
                                                        (
                                                            etalist[0]*Deltaphi
                                                            + betaval*rho1[ind]*(phi[1][ind] - phi[0][ind])
                                                            - rho0[ind] / tauval
                                                            - nt * phi[0][ind]
                                                        );
                    }else if(n==0){
                        aval = 1.0/(lambda_rho + 1.0/tauval) * 
                                                        (
                                                            Dtphi + etalist[0]*Deltaphi
                                                            + betaval*rho1[ind]*(phi[1][ind] - phi[0][ind])
                                                            - rho0[ind] / tauval
                                                            + nt * psi[0][i*n1+j]
                                                        );
                    }else{
                        aval = 1.0/(lambda_rho + 1.0/tauval) * 
                                                        (
                                                            Dtphi + etalist[0]*Deltaphi
                                                            + betaval*rho1[ind]*(phi[1][ind] - phi[0][ind])
                                                            - rho0[ind] / tauval
                                                        );
                    }
                    newrhovalue=cubic_solve(aval, 0, cval);
                    rho0[ind] = fmin(1,fmax(0,newrhovalue));
                    
                }
            }
            // fftps2d->solve_heat_equation(&rho0[n*n1*n2], etalist[0]);
		}
    }

    void update_rho1(const double* rho0,double* rho1, const double* rho2,const double* mx,const double* my,const double* obstacle){

        memcpy(&rho1[0], rhont0tmps[1], n1*n2*sizeof(double));

        double gammaval = gamma;
        double betaval  = beta;
        double tauval = tau[1];

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;
                    
                    double mvalue=0;
                    double Dtphi =0;
                    double Deltaphi=0;

                    calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[1]);
                    double convval3 = rho0[ind] * (phi[1][ind] - phi[0][ind]);
                    double convval4 = (phi[2][ind] - phi[1][ind]);

                    double newrhovalue = 0;
                    double aval=0,cval=0;

                    if(n==nt-1){
                        aval = 1.0/(lambda_rho + 1.0/tauval) * 
                                ( 
                                    - nt * phi[1][ind]
                                    + etalist[1]*Deltaphi
                                    + betaval * convval3 + gammaval * convval4
                                    - rho1[ind] / tauval
                                );                        
                    }else if(n==0){
                        aval = 1.0/(lambda_rho + 1.0/tauval) * 
                            ( 
                                Dtphi + etalist[1]*Deltaphi
                                + betaval * convval3 + gammaval * convval4
                                - rho1[ind] / tauval
                                + nt * psi[1][i*n1+j]
                            );
                    }else{
                        aval = 1.0/(lambda_rho + 1.0/tauval) * 
                                ( 
                                    Dtphi + etalist[1]*Deltaphi
                                    + betaval * convval3 + gammaval * convval4
                                    - rho1[ind] / tauval
                                );
                    }
                    cval = -0.5*alphalist[1]*mvalue*mvalue / (lambda_rho + 1.0/tauval);

                    newrhovalue=cubic_solve(aval, 0, cval);                    
                    rho1[ind]=fmin(1,fmax(0,newrhovalue));    
                }
            }
        }
	}

    void update_rho2(const double* rho0, const double* rho1, double* rho2,const double* mx,const double* my,const double* obstacle){

        memcpy(&rho2[0], rhont0tmps[2], n1*n2*sizeof(double));

        double tauval = tau[2];

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    int ind = n*n1*n2+i*n1+j;

                    if(obstacle[i*n1+j] > 0){
                        rho2[ind]=0;
                    }else{
                        double mvalue=0;
                        double Dtphi =0;
                        double Deltaphi=0;

                        calculate_rho_related(mvalue, Dtphi, Deltaphi, n, i, j, mx, my, phi[2]);

                        double cval = -0.5*alphalist[2]*mvalue*mvalue / (lambda_rho + 1.0/tauval);

                        double aval = 0;
                        if(n==nt-1){
                            aval = (- nt * phi[2][ind] + etalist[2]*Deltaphi - rho2[ind]/tauval) / (lambda_rho + 1.0/tauval);                       
                        }else if(n==0){
                            aval = (Dtphi              + etalist[2]*Deltaphi - rho2[ind]/tauval + nt * psi[2][i*n1+j]) / (lambda_rho + 1.0/tauval);
                        }else{
                            aval = (Dtphi              + etalist[2]*Deltaphi - rho2[ind]/tauval) / (lambda_rho + 1.0/tauval);
                        }

                        double newrhovalue=cubic_solve(aval, 0, cval);
                        rho2[ind]=fmin(1,fmax(0,newrhovalue));
                    }
                }
            }
        }
    }

    double calculate_grad_mx(const double* mxTmp, const int n, const int im, const int i, const int ip, const int jm, const int j, const int jp){
        return n1*(mxTmp[n*n1*n2+i*n1+j]-mxTmp[n*n1*n2+i*n1+jm]);
    }

    double calculate_grad_my(const double* myTmp, const int n, const int im, const int i, const int ip, const int jm, const int j, const int jp){
        return n2*(myTmp[n*n1*n2+i*n1+j]-myTmp[n*n1*n2+im*n1+j]);
    }

    double calculate_dtrho(const double* rho, const int n, const int i, const int j){
    	double dtrho=0;
    	if(n==0) dtrho=0;
        else     dtrho=1.0*nt*(rho[(n)*n1*n2+i*n1+j]-rho[(n-1)*n1*n2+i*n1+j]); 
        return dtrho;
    }

    double calculate_Delta_value(const double* rho, const int n, const int im, const int i, const int ip, const int jm, const int j, const int jp){
        return - n1*n1 * (-rho[n*n1*n2+i*n1+jm]+2*rho[n*n1*n2+i*n1+j]-rho[n*n1*n2+i*n1+jp])
               - n2*n2 * (-rho[n*n1*n2+im*n1+j]+2*rho[n*n1*n2+i*n1+j]-rho[n*n1*n2+ip*n1+j]);
    }

    double update_phi_all(double* const rho[], double* const mx[], double* const my[], const double* obstacle){

    double dtrho0, nablamx0, nablamy0, dtrho1, nablamx1, nablamy1, dtrho2, nablamx2, nablamy2, Deltarho0, Deltarho1, Deltarho2, convval, convval_gamma;
    double dtrho3, nablamx3, nablamy3, Deltarho3;
    int n,ip,im,jp,jm,ind;

    double error = 0;

        for(int n=0;n<nt;++n){  
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){

                    ind = n*n1*n2+i*n1+j;

                    setup_indices(im,ip,jm,jp,i,j);

                    if(n == 0){

                    }else{
                        dtrho0 = calculate_dtrho(rho[0], n, i, j);
                        nablamx0=calculate_grad_mx(mx[0],n,im,i,ip,jm,j,jp);
                        nablamy0=calculate_grad_my(my[0],n,im,i,ip,jm,j,jp);
                        Deltarho0 = calculate_Delta_value(rho[0],n,im,i,ip,jm,j,jp);

                        dtrho1 = calculate_dtrho(rho[1], n, i, j);
                        nablamx1=calculate_grad_mx(mx[1],n,im,i,ip,jm,j,jp);
                        nablamy1=calculate_grad_my(my[1],n,im,i,ip,jm,j,jp);
                        Deltarho1 = calculate_Delta_value(rho[1],n,im,i,ip,jm,j,jp);

                        dtrho2 = calculate_dtrho(rho[2], n, i, j);
                        nablamx2=calculate_grad_mx(mx[2],n,im,i,ip,jm,j,jp);
                        nablamy2=calculate_grad_my(my[2],n,im,i,ip,jm,j,jp);
                        Deltarho2 = calculate_Delta_value(rho[2],n,im,i,ip,jm,j,jp);

                        fftps[0]->u[ind]  = - (dtrho0 + nablamx0 + nablamy0 + beta*rho[0][ind]*rho[1][ind]                     - etalist[0]*Deltarho0 ); 
                        fftps[1]->u[ind]  = - (dtrho1 + nablamx1 + nablamy1 - beta*rho[0][ind]*rho[1][ind] + gamma*rho[1][ind] - etalist[1]*Deltarho1 );
                        fftps[2]->u[ind]  = - (dtrho2 + nablamx2 + nablamy2 - gamma*rho[1][ind]                                - etalist[2]*Deltarho2 ); 

                        
                    }

                    
                    
                }
            }
        }

        fftps[0]->perform_inverse_laplacian(2*beta, etalist[0]);
        fftps[1]->perform_inverse_laplacian(2*beta + 2*gamma, etalist[1]);
        fftps[2]->perform_inverse_laplacian(0, etalist[2]);

        for(int k=0;k<3;++k){
            for(int i=0;i<n1*n2*nt;++i){
                phi[k][i] += sigma[k] * fftps[k]->workspace[i];
                error += fftps[k]->u[i] * fftps[k]->workspace[i];
            }
        }
        return error/(1.0*n1*n2*nt);

    }

    void update_psi(double* const rho[]){
        for(int k=0;k<3;++k){
            for(int i=0;i<n1*n2;++i){
                psi[k][i] += sigma[0] * (rho[k][i] - rhont0tmps[k][i]) * nt;
            }
        }
    }

    double calculate_energy_num(double* const rho[],double* const  mx[], double* const my[], const int num) const{
        double sum1=0;

        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int ind = n*n1*n2+i*n1+j;

                    double mxval=0.5*mx[num][ind];
                    double myval=0.5*my[num][ind];

                    if(j>0) mxval += 0.5*mx[num][n*n1*n2+i*n1+j-1];
                    if(i>0) myval += 0.5*my[num][n*n1*n2+(i-1)*n1+j];

                    double mval=sqrt(mxval*mxval+myval*myval);

                    double rhoval=rho[num][ind];

                    if(rhoval>0){
                        sum1 += alphalist[num]*mval*mval/(2.0*rhoval);
                    }
                }
            }
        }

        return sum1/(n1*n2*nt);
    }


    double calculate_energy(double* const rho[], const double* obstacle) const{
        return -1;
    }

    double calculate_dual(double* const rho[], const double* obstacle) const{   

        double term0=0;
        double term1=0;
        double term2=0;
        double term3=0;
        double term4=0;

        for(int i=0;i<n2;++i){
            for(int j=0;j<n1;++j){
                term0 += phi[0][i*n1+j]*rho[0][i*n1+j] - phi[0][(nt-1)*n1*n2+i*n1+j]*rho[0][(nt-1)*n1*n2+i*n1+j];
                term1 += phi[1][i*n1+j]*rho[1][i*n1+j] - phi[1][(nt-1)*n1*n2+i*n1+j]*rho[1][(nt-1)*n1*n2+i*n1+j];
                term2 += phi[2][i*n1+j]*rho[2][i*n1+j] - phi[2][(nt-1)*n1*n2+i*n1+j]*rho[2][(nt-1)*n1*n2+i*n1+j];
            }
            
        }
        for(int n=0;n<nt;++n){
            for(int i=0;i<n2;++i){
                for(int j=0;j<n1;++j){
                    int idx = n*n1*n2+i*n1+j;
                    term4 += beta   * rho[0][idx] * rho[1][idx] * (phi[0][idx] - phi[1][idx]);
                    term4 += gamma  * rho[1][idx] * rho[2][idx] * (phi[1][idx] - phi[2][idx]);
                }
            }
        }

        return (term0+term1+term2+term3)/(n1*n2) + term4/(n1*n2*nt);
    }

    void update_step_sizes(double& tau, double& sigma, const double error, const double rel_dual_error, const double beta_1, const double beta_2){
    	if(error > rel_dual_error*beta_1){
            tau   *= beta_2;
            sigma /= beta_2;
        }

        if(error < rel_dual_error/beta_1){
            tau   /= beta_2;
            sigma *= beta_2;
        }
    }


    void display_log(const int iterPDHG, const double tau, const double sigma, const double rel_dual_error) const{
        printf("iter: %5d tau: %5.2f sigma: %5.2f rel error: %10.2e\n", iterPDHG+1, tau, sigma, rel_dual_error);
    }

    void run(double* rho[], const double* obstacle, int skip=1){
        for(int k=0;k<3;++k){
            memcpy(rhont0tmps[k], &rho[k][0], n1*n2*sizeof(double));
        }
        double error=1, dual_gap=1, energy=1, dual=0, sanity_value=1, rel_dual_error = 1;
        int iterPDHG;

        double beta_1 = 1.5;
        double beta_2 = 0.9;

        double rel_dual_error_previous = 100;

        for(iterPDHG=0; iterPDHG<max_iteration; ++iterPDHG){
            // update phi
            for(int k=0;k<3;++k){
                memcpy(phitmps[k],phi[k],n1*n2*nt*sizeof(double));
            }
            rel_dual_error  = update_phi_all(rho,mx,my,obstacle);
            rel_dual_error_previous = rel_dual_error;

            for(int k=0;k<3;++k){
                // fftps[0]->solve_heat_equation(rho[k], 1e-5);
                for(int i=0;i<n1*n2*nt;++i){
                    phi[k][i] = 2*phi[k][i] - phitmps[k][i];
                }
            }

            // get the data before updates
            for(int k=0;k<3;++k){
                memcpy(rhotmps[k],rho[k],n1*n2*nt*sizeof(double));
            }

            update_rho0(rho[0],rhotmps[1],rhotmps[2],rhotmps[2],mx[0],my[0],obstacle);
            update_m(mx[0],my[0],rho[0],phi[0],0);    
            update_rho1(rhotmps[0],rho[1],rhotmps[2],mx[1],my[1],obstacle);
            update_m(mx[1],my[1],rho[1],phi[1],1);
            update_rho2(rhotmps[0],rhotmps[1],rho[2],mx[2],my[2],obstacle);
            update_m(mx[2],my[2],rho[2],phi[2],2);

            update_psi(rho);

            if((iterPDHG+1)%skip==0){
                display_log(iterPDHG,  tau[0],  sigma[0],  rel_dual_error);
                create_bin_file(rho[0], n1*n2*nt, "./data/rho0.csv");
                create_bin_file(rho[1], n1*n2*nt, "./data/rho1.csv");
                create_bin_file(rho[2], n1*n2*nt, "./data/rho2.csv");
            }

            // stopping condition
            if((iterPDHG>20 ) && (fabs(rel_dual_error)<tolerance)) break;
        }

        cout<<"The method is done!!"<<endl;
        display_log(iterPDHG,  tau[0],  sigma[0],  rel_dual_error);
    }




}; // Method class


#endif