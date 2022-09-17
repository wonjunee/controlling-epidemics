#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cstring>
// #include "poisson_solver.h"
#include "poisson_solver_3d.h"
#include "helper.h"
#include "initializer.h"
#include "method.h"

using namespace std;

int main(int argc, char **argv)
{

    if(argc!=11){
        cout << "Need to do the following : " << endl;
        cout << argv[0] <<  " [n1] [n2] [nt] [tau] [sigma] [tolerance] [max iteration] [skip] [beta] [gamma]" << endl;
        exit(1);
    }

    // Parameters for Grids
    
    int n1      = stoi(argv[1]);  // x panels
    int n2      = stoi(argv[2]);  // y panels
    int nt      = stoi(argv[3]);  // t panels
    double tau  =stod(argv[4]);
    double sigma=stod(argv[5]);
    double tolerance =stod(argv[6]);
    int max_iteration=stoi(argv[7]);
    int skip=stoi(argv[8]);
    
    double beta  = 0.8;    // infection rate
    double gamma = 0.3;    // recovery rate

    double base=0;

    // coefficients for velocities
    double alphalist[3] = {1, 100, 10000};
    
    // initializing obstacle C
    // 1 if x in C, 0 otherwise.
    double* obstacle = new double[n1*n2];
    for(int i=0;i<n1*n2;++i) obstacle[i] = 0;

    double eta   = 0.0025; // viscosity constants in PDE
    double etalist[] = {eta, eta, eta};  

    // Initialize rho0, rho1, rho2
    double* rho[3];
    for(int i=0;i<3;++i) {
        rho[i] = new double[n1*n2*nt];
    }

    intialize_rho0(rho[0], n1, n2, nt);
    intialize_rho1(rho[1], n1, n2, nt);
    intialize_rho2(rho[2], n1, n2, nt);

    // initialize the method
    Method method(n1, n2, nt, tau, sigma, max_iteration, tolerance, alphalist, etalist, beta, gamma);

    cout<<"XXX G-Prox PDHG XXX"<<endl;
    cout<<endl;
    cout<<"n1 : "<<n1<<" n2 : "<<n2<<" nt : "<<nt<<" base : "<<base<<endl;
    cout<<fixed;
    cout<<"tau           : "<<tau<<endl;
    cout<<"sigma         : "<<sigma<<endl;
    cout<<"max_iteration : "<<max_iteration<<endl;
    cout<<"tolerance     : "<<scientific<<tolerance<<endl;
    cout<<"eta           : "<<scientific<<eta<<endl;
    cout<<"beta          : "<<scientific<<method.beta<<endl;
    cout<<"gamma         : "<<scientific<<method.gamma<<endl;

    create_csv_file_for_parameters(n1,n2,nt,0,0,0,method.beta,method.gamma);

    cout<<"\nXXX Starting Iterations XXX"<<endl;

    clock_t t;
    t = clock();

    method.run(rho,obstacle,skip);   

    t = clock() - t;

    printf ("CPU time for Iterations: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    for(int k=0;k<3;++k){
        delete[] rho[k];
    }
    delete[] obstacle;
}