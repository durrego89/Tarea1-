#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define MAX(x,y) x>y?x:y

void w2u(double *w, double gamma, double *u);
void u2w(double *u, double gamma, double *w);

void euler_flux(double *w, double gamma, double *f);

void roe_flux(double *wL, double *wR, double gamma);

double calc_time_step(double cfl,double  dx,double  gamma,int  bcells, double  *U);
void update_solution(double **U, double *fluxes, double dt, double dx, double gamma, int bcells);



int main() {
   

    double gamma=1.4;
    int cells = 1000;
    int bcells = cells+2;
    double dx =1.0/cells;
    
    double cfl = 0.8;
    double t = 0.0;
    double tf = 0.2;
    double nsteps = 0;
    
    double pasos = 50000;
    
    int i;
    int N = 5000;
    double **U = NULL;
    double **fluxes = NULL;
   
    

    //void allocate(){
        int j;
        
        U = malloc(N * sizeof(double *));
        fluxes = malloc(N * sizeof(double *));

        for(j = 0; j < N; j++){
            U[j] = malloc(3 * sizeof(double));
            fluxes[j] = malloc(3 * sizeof(double));

        }
    //}
    
    test_case_sod(U, bcells, gamma);
//    
//    for(i=1;i<=pasos;i++){
//        if(t==tf){
//            break;
//        }
//        dt = calc_time_step(cfl, dx, gamma, bcells, U);
//        if (t+dt > tf){
//            dt = tf - t;
//        }
//        update_solution(U, fluxes, dt, dx, gamma, bcells);
//        t += dt;
//        nsteps += 1;
// 
//    }

    
    
    return 0;
}
//==========================================================================================

void test_case_sod(double **U, double bcells, double gamma){

//    Populate the initial data with the Sod test case.

    int mid_point = int(bcells/10.0*3.0);
    double sod_L[2];
    sod_L[0]=1.0;
    sod_L[1]=0.75;
    sod_L[2]=1.0;
    
    double sod_R[2];
    sod_R[0]=0.125;
    sod_R[1]=0.0;
    sod_R[2]=0.10;

    int i;

//    for(i=0,i<=mid_point;i++){
//        U[i][j] = w2u(sod_L, gamma);  //REVISAR
//
//        U[i][j] = w2u(sod_R, gamma);
    //}
    }
//==========================================================================================

//void update_solution(double *U, double *fluxes, double dt, double dx, double gamma, int bcells){
//
////Updates the solution of the equation
////via the Godunov procedure.
//    
//    double *f1;
//      double *f2;
//    double n_x=3;
//    double wL, wR;
//    f1 = malloc (n_x * sizeof(double));
//    f2 = malloc (n_x * sizeof(double));
//    
//    int i;
//    
//// Create fluxes
//    for(i=1;i<=bcells;i++){
//        u2w(U, gamma,f1); // REVISAR u[i]
//        wL =f1[i];
//        u2w(U, gamma,f2); //REVISAR Ui+1
//        wR =f2[i];
//        //fluxes[i] = roe_flux(wL, wR, gamma); // REVISAR
//    }
//// Update solution
//    for(i=1;i<=bcells;i++){
//        U[i] +=  (dt/dx) * (fluxes[i-1]-fluxes[i]);
//    }
//                            
//    U[0] = U[1];
//    U[bcells-1] = U[bcells-2];
//}
//
////===========================================================================================
//
//double calc_time_step(double cfl,double  dx,double  gamma,int  bcells,double  *U){
//
////Calculates the maximum wavespeeds and thus the timestep
////via an enforced CFL condition.
//
//    double max_speed = -1.0;
//    int i;
//    double  c, u;
//    double w[bcells-1];
//    double dt;
//    double *f;
//    double n_x=3;
//    
//    
//    f = malloc (n_x * sizeof(double));
//   
//    for(i=1;i<bcells-1;i++){
//        u2w(U, gamma,f);  // REVISAR PASO EXTRANO
//        w[i] =f[i] ;
//        u = w[1];
//        c = sqrt(gamma*w[2]/w[0]);
//        max_speed = MAX(max_speed, abs(u)+c);
//    }
//    dt = cfl*dx/max_speed;  // CFL condition;
//    return dt;
//    }
//
//
////==============================================================================================================
//
//void roe_flux(double *wL, double *wR, double gamma){
//
//    double *uL, *uR;
//    double *f1;
//    double *f2;
//    double *u;
//    double *flux;
//    
//    double n_x=3;
//    
//    uR = malloc (n_x * sizeof(double));
//    uL = malloc (n_x * sizeof(double));
//    flux = malloc (n_x * sizeof(double));
//    
//    f1 = malloc (n_x * sizeof(double)); // auxiliar
//    f2 = malloc (n_x * sizeof(double)); // auxiliar
//    u = malloc (n_x * sizeof(double));  // auxiliar
//    
//    double rho, rhoL, rhoR;
//    double v, vL, vR;
//    double pL, pR;
//    double a, aL, aR;
//    double H, HL, HR;
//    double RT;
//    
//    
//    int j,i;
//    
////Use the Roe approximate Riemann solver to calculate fluxes.
//
// 
//    
//    
//     w2u(wL, gamma ,u);
//    for(j=0;j<n_x;j++){
//       uL[j]=u[j];
//    }
//    
//    w2u(wR, gamma ,u);
//    for(j=0;j<n_x;j++){
//        uL[j]=u[j];
//    }
//
//// Primitive and other variables.
//// Left state
//    rhoL = wL[0];
//    vL = wL[1];
//    pL = wL[2];
//    aL = sqrt(gamma*pL/rhoL); //revisar
////    HL = ( uL[2] + pL ) / rhoL;
//
//// Right state
//    rhoR = wR[0];
//    vR = wR[1];
//    pR = wR[2];
//    aR = sqrt(gamma*pR/rhoR);
//    HR = ( uR[2] + pR ) / rhoR;
//
//// First compute the Roe Averages
//    RT = sqrt(rhoR/rhoL);
//    rho = RT*rhoL;
//    v = (vL+RT*vR)/(1.0+RT);
//    H = (HL+RT*HR)/(1.0+RT);
//    a = sqrt((gamma-1.0)*(H-0.5*v*v));
//
//// Differences in primitive variables.
//    double drho, du, dP;
//    
//    drho = rhoR - rhoL;
//    du = vR - vL;
//    dP = pR - pL;
//
//// Wave strength (Characteristic Variables).
//    double *dV;
//    dV = malloc (n_x * sizeof(double));
//    
//    dV[0] = 0.5*(dP-rho*a*du)/(a*a);
//    dV[1] = -( dP/(a*a) - drho );
//    dV[2] = 0.5*(dP+rho*a*du)/(a*a);
//
//// Absolute values of the wave speeds (Eigenvalues)
//    double *ws;
//    ws = malloc (n_x * sizeof(double));
//    
//
//    ws[0] = abs(v-a);
//    ws[1] = abs(v);
//    ws[2] = abs(v+a);
//
//// Modified wave speeds for nonlinear fields (the so-called entropy fix, which
//// is often implemented to remove non-physical expansion shocks).
//// There are various ways to implement the entropy fix. This is just one
//// example. Try turn this off. The solution may be more accurate.
//  
//    double Da;
//    Da = MAX(0.0, 4.0*((vR-aR)-(vL-aL)) );
//    if (ws[0] < 0.5*Da){
//        ws[0] = ws[0]*ws[0]/Da + 0.25*Da;
//    }
//    Da = MAX(0.0, 4.0*((vR+aR)-(vL+aL)) );
//    if (ws[2] < 0.5*Da){
//        ws[2] = ws[2]*ws[2]/Da + 0.25*Da;
//    }
//// Right eigenvectors
//    double R[3][3];
//    
//    
//    R[0][0] = 1.0;
//    R[1][0] = v - a;
//    R[2][0] = H - v*a;
//
//    R[0][1] = 1.0;
//    R[1][1] = v;
//    R[2][1] = 0.5*v*v;
//
//    R[0][2] = 1.0;
//    R[1][2] = v + a;
//    R[2][2] = H + v*a;
//
//// Compute the average flux.
//    
//    euler_flux(wL, gamma,f1);
//    euler_flux(wR, gamma,f2);
//    
//    for(i=0;i<n_x;i++){
//    flux[i] = 0.5*(f1[i]  + f2[i]  );
//    }
//    
//// Add the matrix dissipation term to complete the Roe flux.
//    
//    
//
//    for(i=0;i<n_x;i++){
//        for(j=0;j<n_x;j++){
//
//            flux[i] -= 0.5*ws[j]*dV[j]*R[i][j];
//
//        }
//    }
//}
//
////==============================================================================================================
//
//
//void euler_flux(double *w, double gamma, double *f){
//
////Calculate the conservative Euler fluxes.
//
//    double rho;
//    double u;
//    double p;
//    double a2;
//    
//    double f_1, f_2, f_3;
//    
//    rho = w[0];
//    u = w[1];
//    p = w[2];
//
//    a2 = gamma*p/rho;
//
//    f[0] = rho*u;  //f_1
//    f[1] = rho*u*u + p;  //f_2
//    f[2] = rho*u*( a2/(gamma-1.0) + 0.5*u*u ); //f_3
//}
//
//
//
//
//
//
////----------------------------------------------------------------------------------------------------------
//
//
//void w2u(double *w, double gamma, double *u){
//    u[0] = w[0];
//    u[1] = w[0]*w[1];
//    u[2] = w[2]/(gamma-1.0)+0.5*w[0]*w[1]*w[1];
//}
//
//
//
//void u2w(double *u, double gamma, double *w){
//
//    w[0] = u[0];
//    w[1] = u[1]/u[0];
//    w[2] = (gamma-1.0)*( u[2] - 0.5*w[0]*w[1]*w[1] );
//}
//
//
//
//
//
//
//



















