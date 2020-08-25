#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//Declare pi
#define M_PI 3.14159265358979323846
/*Delta N_HI*/
#define DNHI 2.5e16
/*N_grid*/
#define NGRID 2000
/*Delta t^prime == 2.5e15 cm^-2*/
#define DTIMESTEP 2.5e15
/*N t^prime evolve until t'=3*10^19 cm^-2*/
#define NTIMESTEP 12000

/*Number of Frequency bins*/
#define N_NU 128
/*Rydberg temperature*/
#define RYD_K 157800
/*Helium ionization energy*/
#define ION_HE 1.809
/*f_He*/
#define ABUND_HE 0.079

//define the demension of energy transfering matrix
#define N 5

/* Sets blackbody incident spectrum at given T */
void set_bb(double *frac, double T) {
  int i;
  double nu, tot=0.;
  double myxs = 0, epsilon = 0;

  T /= RYD_K;
  for(i=0;i<N_NU;i++) {
    nu = pow(4., (i+0.5)/N_NU);
#ifdef THIS_HCOL
    epsilon = sqrt(nu-1.);
    myxs = THIS_HCOL * 6.30e-18*pow(nu,-4.)*exp(4.-4.*atan2(epsilon,1.)/epsilon)/(1.-exp(-2*M_PI/epsilon));
#endif
    tot += frac[i] = nu*nu*nu/(exp(nu/T)-1.) * exp(-myxs);
  }
  for(i=0;i<N_NU;i++) frac[i]/=tot;
}

/* Sets cross sections [in cm^2] */
void set_sigma(double *nu, double *sigH, double *sigHe) {
  int i;
  double epsilon, x, y;

  for(i=0;i<N_NU;i++) {
    nu[i] = pow(4., (i+0.5)/N_NU);
    epsilon = sqrt(nu[i]-1.);
    sigH[i] = 6.30e-18*pow(nu[i],-4.)*exp(4.-4.*atan2(epsilon,1.)/epsilon)/(1.-exp(-2*M_PI/epsilon));
    sigHe[i] = 0.;
    if (nu[i]>ION_HE) {

      /* Verner et al. 1996
       * sigma = sigma_0 [ (x-1)^2 + y_w^2 ] y^{0.5P-5.5} ( 1 + sqrt{y/y_a} )^{-P}
       * y = sqrt(x^2 + y_1^2)
       * x = E/E_0-y_0
       */
      x = nu[i]/1.0007-0.4434;
      y = sqrt(x*x+2.136*2.136);
      sigHe[i] = 9.492e-16 * ( (x-1.)*(x-1.)+2.039*2.039 ) * pow(y,0.5*3.188-5.5) * pow(1.+sqrt(y/1.469), -3.188);
    }
#if 0
    sigH[i] = 6.2e-18*pow(nu[i],-3.);
    sigHe[i] = 0.;
    if (nu[i]>ION_HE) sigHe[i] = 7.9e-18*pow(nu[i]/ION_HE,-3.);
#endif
  }
}

/* Takes in neutral fractions y1H[] and y1He[], and the incident photon spectrum.
 * Outputs rates: delta(y1H), delta(y1He), delta(EH) [Ryd/H]
 *  are done with a "time variable" that is 1 photon/cm^2
 *  [i.e. scaled time = F * t', where F = ionizing photon flux in ph/cm^2/s]
 */
void get_ion_rate(double *y1H, double *y1He, double *fracflux, double *dy1H, double *dy1He, double *dEH, double *Te, int istep, double U) {
  static int is_initialized = 0;
  static double *nu, *sigH, *sigHe;
  long i,j;
  double tauH, tauHe, flux, tautot, wt;

  /* Build internal tables */
  if (!is_initialized) {
    is_initialized = 1;
    /*Memory allocation, sigH/sigHe/nu are all address*/
    sigH = (double*)malloc((size_t)(3*N_NU*sizeof(double)));
    sigHe = sigH + N_NU;
    nu = sigHe + N_NU;
    /*Why set_sigma again here??*/
    set_sigma(nu,sigH,sigHe);
  }

  /*Set initial condition for the change of aboundance and energy*/
  for(j=0;j<NGRID;j++) {
    dy1H[j] = dy1He[j] = dEH[j] = 0.;
  }
  /* Now get contribution from each bin */
  for(i=0;i<N_NU;i++) {
    flux = fracflux[i];
    for(j=0;j<NGRID;j++) {
      /* Optical depth in this slice */
      tauH = DNHI * sigH[i] * y1H[j] / (1 - U/3.e10);
      tauHe = ABUND_HE * DNHI * sigHe[i] * y1He[j]/ (1 - U/3.e10);
      tautot = tauH + tauHe;
//      tautot = tautot / cos(M_PI/4.);
      /* weight = (mean flux in this slice) */
      /*IMPORTANT comments by Chenxiao: approximation made here for exponential function,
       * flux is affected by the optical depth*/
      wt = flux * (tautot>1e-5? (1.-exp(-tautot))/tautot: 1.-tautot/2.);

      dy1H[j] -= wt * sigH[i] * y1H[j];
      dy1He[j] -= wt * sigHe[i] * y1He[j];

      dEH[j] += wt * sigH[i] * y1H[j] * (nu[i]-1.);
      dEH[j] += wt * sigHe[i] * y1He[j] * (nu[i]-ION_HE) * ABUND_HE;
      flux *= exp(-tautot);
    }
  }


}

/* Cooling rate [divided by n_H] in Ry cm^3/s at given neutral fractions
 */
double get_cooling_rate(double Te, double y1H, double y1He) {

  double q12, q13, gamma;

  /* 1->2 */
  gamma = 0.531+2.71e-5*Te-3.22e-10*Te*Te+3.88e-16*Te*Te*Te;
  if (Te>2.5e4)
    gamma = 0.637+1.47e-5*Te-5.92e-12*Te*Te-3.78e-18*Te*Te*Te;
  q12 = 8.6287e-6/2./sqrt(Te)*gamma*exp(-0.75*RYD_K/Te);

  /* 1->3 */
  gamma = 0.350-2.62e-7*Te-8.15e-11*Te*Te+6.19e-15*Te*Te*Te;
  if (Te>2.5e4)
    gamma = 0.276+4.99e-6*Te-8.85e-12*Te*Te+7.18e-18*Te*Te*Te;
  q13 = 8.6287e-6/2./sqrt(Te)*gamma*exp(-0.8888888889*RYD_K/Te);

  return( (3./4.*q12 + 8./9.*q13)*y1H*(1.-y1H+ABUND_HE*(1.-y1He)) );
}

//Calculate the cofactor of mat[x][y] in temp[][]
void coFactor(double mat[N][N], double temp[N][N], int p, int q, int n) {
    int i=0, j=0;
    for (int row=0; row<n; row++) {
        for (int col = 0; col<n; col++) {
           if (row!=p && col!=q) {
                temp[i][j++] = mat[row][col];
                if (j==n-1) {
                    j=0;
                    i++;
                }
           }
        }
    }
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of mat[][]. */
double determinant(double mat[N][N], int n) {
    double D=0;
    if (n==1)
        return mat[0][0];
    double temp[N][N];
    int sign = 1;
    for (int f=0;f<n;f++) {
        coFactor(mat, temp, 0, f, n);
        D += sign * mat[0][f] * determinant(temp, n-1);
        sign = -sign;
    }
    return D;
}

void inverseMat(double mat[N][N], double temp[N][N]) {
    int i=0, j=0;
    double det = determinant(mat, N);
    double blanck[N][N];
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            coFactor(mat, blanck, j, i, N);
            temp[i][j] = pow(-1, i+j) * determinant(blanck, N-1);
            temp[i][j] /= det;
        }
    }
}

/* function for displaying the matrix */
void display(double mat[N][N], int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
            printf("  %lf", mat[i][j]);
        printf("\n");
    }
}

int main(int argc, char **argv) {
  //Iteration index
  int i;
  long j, istep;
  //Define flux in each bin and cross section
  double fracflux[N_NU], nu[N_NU], sigH[N_NU], sigHe[N_NU];
  //Define H, He neutral fractions
  double y1H[NGRID], y1He[NGRID];
  double dy1H[NGRID], dy1He[NGRID];

  double dEH[NGRID];
  //Define energy and temperatures for five species
  double EH[NGRID],EHI[NGRID],EHII[NGRID],EHeI[NGRID],EHeII[NGRID];
  double Te[NGRID],THI[NGRID],THII[NGRID],THeI[NGRID],THeII[NGRID];
  //Define energy transfering matrix
  double M[N][N];
  //Define energy transferring rate
  double nuTeTHII, nuTeTHI, nuTHITe, nuTHIITHI, nuTeTHeI, nuTHeITe, nuTHIITHeI, nuTHITHeI, nuTeTHeII, nuTHIITHeII, nuTHITHeII, nuTHeITHeII;
  //Define the inverse of M
  double I[N][N];
  //Define blackbody incident temperature and ionization front velocity
  double T, U;

  /*The function sscanf returns an integer which is equal to the number of parameters that were successfully converted*/
  sscanf(argv[1], "%lf", &T);
  set_bb(fracflux,T);
  set_sigma(nu,sigH,sigHe);
#if 0
  for(i=0;i<N_NU;i++)
    printf("%2d %8.6lf %8.6lf %11.5lE %11.5lE\n", i, nu[i], fracflux[i], sigH[i], sigHe[i]);
#endif
  i=0; /* to avoid complaining */

  /* Set up initial conditions */
  for(j=0; j<NGRID; j++) {
    y1H[j] = y1He[j] = 1.-1.e-10;
    EH[j] = EHI[j] =  EHeI[j] = 1.e-30;
    EHII[j] = 1.e-30;
    EHeII[j] = 1.e-30;
  }


  /*U is the ionization front speed???*/
  sscanf(argv[2], "%lf", &U);
  printf("InverseMatrix00,   InverseMatrix01,   InverseMatrix02,   InverseMatrix03,   InverseMatrix04\n");
  for(istep=0;istep<NTIMESTEP;istep++) {
    if (istep==0) {
      for(j=0;j<NGRID;j++) {
         /*New temperatures here*/
        Te[j] = EH[j]/1.5/(1.-y1H[j]+ABUND_HE*(1.- 1*y1He[j]))*RYD_K;
        THI[j] = EHI[j]/1.5/y1H[j]*RYD_K;
        THII[j] = EHII[j]/1.5/(1.-y1H[j])*RYD_K;
        THeI[j] = EHeI[j]/1.5/(ABUND_HE*y1He[j])*RYD_K;
        THeII[j] = EHeII[j]/1.5/(ABUND_HE*(1.-y1He[j]))*RYD_K;
      }
     }
    get_ion_rate(y1H,y1He,fracflux,dy1H,dy1He,dEH,Te,istep,U);
   
   //Ub stands for relativistic correction
    double Ub = U / (1 - U/3.e10);
    
    for(j=0;j<NGRID;j++)
      dEH[j] -= get_cooling_rate(Te[j], y1H[j], y1He[j])/(1.+ABUND_HE)/Ub;
    for(j=0;j<NGRID;j++) {
      y1H[j] += DTIMESTEP * dy1H[j];
      y1He[j] += DTIMESTEP * dy1He[j];

      EH[j] += DTIMESTEP * dEH[j];
      Te[j] = EH[j]/1.5/(1.-y1H[j]+ABUND_HE*(1.- 1*y1He[j]))*RYD_K;

      if (istep >= 0) {
      //Set up energy transferring matrix in the column order of: EH, EHII


    nuTeTHII = 0.0753764/Ub*pow(Te[j],-3./2)*1.4;
    
    if (nuTeTHII > 1.)
        nuTeTHII = 1.;
    nuTeTHI = 1/Ub*393.542*(1.08779e-14 - 3.51267e-15*pow(Te[j],0.0902082))*pow(Te[j],1./2);
    nuTHITe = nuTeTHI;
    nuTHIITHI = 1/Ub*8424.81*(7.99395e-15 + 2.51304e-13*pow(THII[j],-0.634768))*pow(THII[j],1./2);

    nuTeTHeI = 1/Ub*98.9781*(6.87675e-16 - 2.18539e-21*pow(Te[j],0.984171))*pow(Te[j],1./2);
    nuTHeITe = nuTeTHeI;
    nuTHIITHeI = 1.31603e-9/Ub;
    nuTHITHeI = 1/Ub*8424.81*(-3.53379e-14 + 3.89354e-14*pow(THI[j],-0.0086412))*pow(THI[j],1./2);

    nuTeTHeII = 0.0758305/Ub*pow(Te[j],-3./2)*1.4;
    if (nuTeTHeII > 1.)
        nuTeTHeII = 1.;
    nuTHIITHeII = 1.40535/Ub*pow(0.143916*THeII[j] + 0.572216*THII[j],-3./2)*1.4;
    if (nuTHIITHeII > 1.)
        nuTHIITHeII = 1.;
    nuTHITHeII = 2.3739e-9/Ub;
    nuTHeITHeII = 1/Ub*4225.07*(7.71629e-11-7.71569e-11*pow(THeII[j],4.18115e-6))*pow(THeII[j],1./2);

     //Use THII because I derive nuTHIIHeII from nuTHIIHI

     //g factor is equivalent to f(\alpha) in the paper, to distinguish from fHe.

     double ge = (1.-y1H[j]+ABUND_HE*(1.- 1*y1He[j]));
     double gHII = (1.-y1H[j]);
     double gHI = y1H[j];
     double gHeI = (ABUND_HE*y1He[j]);
     double gHeII = (ABUND_HE*(1.-y1He[j]));

     M[0][0] = 1. + (nuTeTHII * gHII + nuTeTHI * gHI + nuTeTHeI * gHeI + nuTeTHeII * gHeII) * DTIMESTEP;
     M[0][1] = - nuTeTHII * DTIMESTEP * gHII;
     M[0][2] = - nuTeTHI * DTIMESTEP * gHI;
     M[0][3] = - nuTeTHeI * DTIMESTEP * gHeI;
     M[0][4] = - nuTeTHeII * DTIMESTEP * gHeII;

     M[1][0] = - nuTeTHII * DTIMESTEP * ge;
     M[1][1] = 1. + (nuTeTHII * ge + nuTHIITHI * gHI + nuTHIITHeI * gHeI + nuTHIITHeII * gHeII) * DTIMESTEP;
     M[1][2] = - nuTHIITHI * DTIMESTEP * gHI;
     M[1][3] = - nuTHIITHeI * DTIMESTEP * gHeI;
     M[1][4] = - nuTHIITHeII * DTIMESTEP * gHeII;

     M[2][0] = - nuTHITe * DTIMESTEP * ge;
     M[2][1] = - nuTHIITHI * DTIMESTEP * gHII;
     M[2][2] = 1. + (nuTHITe * ge + nuTHIITHI * gHII + nuTHITHeI * gHeI + nuTHITHeII * gHeII) * DTIMESTEP;
     M[2][3] = - nuTHITHeI * DTIMESTEP * gHeI;
     M[2][4] = - nuTHITHeII * DTIMESTEP * gHeII;

     M[3][0] = - nuTHeITe * DTIMESTEP * ge;
     M[3][1] = - nuTHIITHeI * DTIMESTEP * gHII;
     M[3][2] = - nuTHITHeI * DTIMESTEP * gHI;
     M[3][3] = 1. + (nuTHeITe * ge + nuTHIITHeI * gHII + nuTHITHeI * gHI + nuTHeITHeII * gHeII) * DTIMESTEP;
     M[3][4] = - nuTHeITHeII * DTIMESTEP * gHeII;

     M[4][0] = - nuTeTHeII * DTIMESTEP * ge;
     M[4][1] = - nuTHIITHeII * DTIMESTEP * gHII;
     M[4][2] = - nuTHITHeII * DTIMESTEP * gHI;
     M[4][3] = - nuTHeITHeII * DTIMESTEP * gHeI;
     M[4][4] = 1. + (nuTeTHeII * ge + nuTHIITHeII * gHII + nuTHITHeII * gHI + nuTHeITHeII * gHeI) * DTIMESTEP;

     //Calculate the inverse of M, here is the identity matrix minus the energy transfering matrix
     inverseMat(M, I);

     if (j==13){
/*        
        printf("step=%ld ", istep);
        printf("M13\n");
        display(M, 5, 5);
        printf("I13\n");
        display(I, 5, 5);
        printf("nuTeTHII = %12.5e\n", nuTeTHII);
        printf("nuTeTHI = %12.5e\n", nuTeTHI);
        printf("nuTHIITHI = %12.5e\n", nuTHIITHI);
        printf("nuTeTHeI = %12.5e\n", nuTeTHeI);
        printf("nuTHIITHeI = %12.5e\n", nuTHIITHeI);
        printf("nuTHITHeI = %12.5e\n", nuTHITHeI);
        printf("nuTeTHeII = %12.5e\n", nuTeTHeII);
        printf("nuTHIITHeII = %12.5e\n", nuTHIITHeII);
        printf("nuTHITHeII = %12.5e\n", nuTHITHeII);
        printf("nuTHeITHeII = %12.5e\n", nuTHeITHeII);
  */     
        printf("Te = %12.5e\n", Te[j]);
    /*    
        printf("THII = %12.5e\n", THII[j]);
        printf("THI = %12.5e\n", THI[j]);
        printf("THeI = %12.5e\n", THeI[j]);
        printf("THeII = %12.5e\n", THeII[j]);
        printf("\n");
      */  
     }
        double Tej = Te[j];
        double THIIj = THII[j];
        double THIj = THI[j];
        double THeIj = THeI[j];
        double THeIIj = THeII[j];
        Te[j]   = I[0][0]*Tej + I[0][1]*THIIj + I[0][2]*THIj + I[0][3]*THeIj + I[0][4]*THeIIj;
        THII[j] = I[1][0]*Tej + I[1][1]*THIIj + I[1][2]*THIj + I[1][3]*THeIj + I[1][4]*THeIIj;
        THI[j]  = I[2][0]*Tej + I[2][1]*THIIj + I[2][2]*THIj + I[2][3]*THeIj + I[2][4]*THeIIj;
        THeI[j] = I[3][0]*Tej + I[3][1]*THIIj + I[3][2]*THIj + I[3][3]*THeIj + I[3][4]*THeIIj;
        THeII[j]= I[4][0]*Tej + I[4][1]*THIIj + I[4][2]*THIj + I[4][3]*THeIj + I[4][4]*THeIIj;

        EH[j] = Te[j]*1.5*(1.-y1H[j]+ABUND_HE*(1.- 1*y1He[j]))/RYD_K;
        EHII[j] = THII[j]*1.5*(1.-y1H[j])/RYD_K;
        EHI[j] = THI[j]*1.5*y1H[j]/RYD_K;
        EHeI[j] = THeI[j]*1.5*ABUND_HE*y1He[j]/RYD_K;
        EHeII[j] = THeII[j]*1.5*ABUND_HE*(1-y1He[j])/RYD_K;
    }
    }
  }

  /*Print out the overall one-dimentional model values for each cell*/
  //n: run the code on cluster, plot dEH[3] and dEH[1245]
  for(i=0; i<5; i++)
    printf("\n");
  printf("One-dimensional model\n");
  printf("Timestep=%7d\n", NTIMESTEP);
  printf("j, (j+.5)*DNHI,     y1H[j],     y1He[j],     EH[j],     Te[j],     EHII[j],     THII[j],     EHI[j],     THI[j],     EHeI[j],     THeI[j],     EHeII[j],     THeII[j],     dEH[j]\n");
  for(j=0; j<NGRID; j++) {
    printf("%4ld %11.5lE %8.15lf %8.6lf %8.6lf %7.15lf %8.6lf %7.15lf %8.6lf %7.15lf %8.6lf %7.15lf %8.6lf %7.15lf %8.10lf\n",
      j, (j+.5)*DNHI, y1H[j], y1He[j], EH[j], Te[j], EHII[j], THII[j], EHI[j], THI[j], EHeI[j], THeI[j], EHeII[j], THeII[j], dEH[j]);
  }

  return(0);
}
