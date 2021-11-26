/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 1001969
 *   Name        : Goh Xing Yang
 *
 ***************************************************************************/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

/* Function to run root finding on shock wave 
*/
void shockwave(const char *q1_file) {
	machdata_t* machDataLine = (machdata_t *)malloc(sizeof(machdata_t));
	assert(machDataLine != NULL);
    FILE *shockFile = safe_fopen(q1_file, "r");
	int numMach = 0;
	double *machArray = (double *)malloc(sizeof(double)*INITIALBUFFER);
	assert(machArray != NULL);

    /* Reading in and solving shock-wave equation */
	read_shockwave_a(shockFile, machDataLine);
	machArray = read_mach_array(shockFile, machArray, &numMach);
    print_shockwave(machDataLine, machArray, numMach);

    /* Free variables and close files */
    free(machArray);
    free(machDataLine);
    fclose(shockFile);
	return;

}
/**************************************************************************/

/* Function to safely open a file
*/
FILE *safe_fopen(const char* path, const char* mode) {
    FILE* fp = fopen(path, mode);
    if (fp == NULL) {
        perror("file open error.");
        exit(EXIT_FAILURE);
    }
    return fp;
}
/**************************************************************************/

/* Function to print an error message and exit the program 
*/
void invalid_data_error(char *msg, int line) {
    printf("Error on line %3d: %s\n", line, msg);
    printf("The program will now terminate\n");
    exit(EXIT_FAILURE);
    return;
}
/**************************************************************************/

/* Checking initial file format and removing the header line then reading
the first section of the data in shockwave file
*/
void read_shockwave_a(FILE *shockFile, machdata_t* machDataLine) {

    char ch;
    double val;

    /* Checking if the data file begins with a number, terminating if so */
    if (fscanf(shockFile, "%lf", &val) == 1) {
        invalid_data_error("Data formatting error", __LINE__);
    }

    /* Taking the header line from the csv file and throwing it away */
    fscanf(shockFile, "%c", &ch);
    while (ch != '\n') {
        fscanf(shockFile, "%c", &ch);
    }

    /* Scanning in the first section of the input file */
    if (fscanf(shockFile, "%lf,%lf,%lf,%lf,%lf\n", &machDataLine->M, 
                            &machDataLine->theta, &machDataLine->beta_u, 
                            &machDataLine->beta_l, &machDataLine->gamma)!=NUMFIELDSSHOCK) {
    	invalid_data_error("Data formatting error", __LINE__);
   }

    /* Removing header for second section */
    fscanf(shockFile, "%c", &ch);
    while (ch != '\n') {
        fscanf(shockFile, "%c", &ch);
    }
   return;
}
/**************************************************************************/

/* Dynamically read the Mach numbers for section b of shockwave file
*/
double *read_mach_array(FILE *shockFile, double *machArray, int *numMach) {
	int buffer = INITIALBUFFER, read = 0;

	/* Scan in the mach numbers one at a time and reallocating when needed */
	while (fscanf(shockFile, "%lf\n", &machArray[read]) == 1) {
        read++;
		if (read == buffer) {
			buffer *= 2;
			machArray = (double *)realloc(machArray, sizeof(double)*buffer);
			assert(machArray != NULL);
		}
	}

    /* Reallocate the mach array to the appropriate size */
    machArray = (double *)realloc(machArray, sizeof(double)*read);
    assert(machArray != NULL);
    *numMach = read;
	return machArray;
}
/**************************************************************************/

/* Function of beta
*/
double f_beta(double M, double theta, double gamma, double beta) {
    return (2*(1/tan(beta))*((M*M*(sin(beta))*(sin(beta))-1)/
           (M*M*(gamma + cos(2*beta))+2))-tan(theta));
}
/**************************************************************************/

/* Derivative of the function of beta 
*/
double f_beta_prime(double M, double theta, double gamma, double beta) {
    return (4*M*M*cos(beta)*sin(beta))/(tan(beta)*
           ((gamma + cos(2*beta))*M*M + 2))
           - (2*(tan(beta)*tan(beta) + 1)*(M*M*sin(beta)*sin(beta) - 1))/
           (tan(beta)*tan(beta)*((gamma + cos(2*beta))*M*M + 2)) +
           (4*M*M*sin(2*beta)*(M*M*sin(beta)*sin(beta) - 1))/
           (tan(beta)*(M*M*(gamma + cos(2*beta)) + 2)*
           (M*M*(gamma + cos(2*beta)) + 2));

}
/**************************************************************************/

/* Newton Raphson root convergence method
*/
double newtonraph(double guess, double theta, double M, double gamma) {
    double fx1, fpx1;
    int iterations=0;
    double beta = guess*M_PI/PI_DEGREES;
    fx1 = f_beta(M, theta*M_PI/PI_DEGREES, gamma, beta); 
    fpx1= f_beta_prime(M, theta*M_PI/PI_DEGREES, gamma, beta);
    
    /* Newton Raphson iteration steps */
    while (fabs(fx1) > EPS && iterations < MAXITER) {
        iterations++;
        beta = beta-(fx1/fpx1);
        fx1 = f_beta(M, theta*M_PI/PI_DEGREES, gamma, beta);
        fpx1= f_beta_prime(M, theta*M_PI/PI_DEGREES, gamma, beta);
    }

    /* Unable to converge to a root */
    if (iterations == MAXITER && fabs(fx1) > EPS) {
        return DETATCHED;
    }
    return beta*PI_DEGREES/M_PI;
}
/**************************************************************************/

/* Printing the shockwave theta-beta-M data to csv file
*/
void print_shockwave(machdata_t* machDataLine, double *machArray, int numMach) {
    FILE *writeFile = safe_fopen("out_shock.csv", "w");
    int i;
    double bl, bu, theta;
    fprintf(writeFile, "M,theta,beta_lower,beta_upper\n");
    for (i = 0; i < numMach; i++) {

        /* Initial guesses for theta = 0 */
        theta = 0;
        bl = asin(1.0/machArray[i])*PI_DEGREES/M_PI;
        bu = THETAMAX;

        /* Calculate angles for increments of theta */
        bl = newtonraph(bl, theta, machArray[i], machDataLine->gamma);
        bu = newtonraph(bu, theta, machArray[i], machDataLine->gamma);

        while (bl!=DETATCHED && bu!=DETATCHED) {
            fprintf(writeFile, "%0.4lf,%0.4lf,%0.4lf,%0.4lf\n", machArray[i], 
                                theta, bl, bu);
            theta++;

            bl = newtonraph(bl, theta, machArray[i], machDataLine->gamma);
            bu = newtonraph(bu, theta, machArray[i], machDataLine->gamma);

            /*  Condition to break the loop */
            if (bl<THETAMIN || bu>THETAMAX || bu<bl) {
                break;
            }

        }
    }
    fclose(writeFile);
    return;
}
/**************************************************************************/

/* Function to solve the tri-diagonal linear system 
*/
void linsolve(const char* q3_file) {
    tridata_t* triArray = (tridata_t *)malloc(sizeof(tridata_t)*INITIALBUFFER);
    assert(triArray != NULL);
    FILE *triSysFile = safe_fopen(q3_file, "r");
    int numRow = 0;

    /* Solving the tri-diagonal linear system using the Thomas algorithm*/
    triArray = read_tri_data(triSysFile, triArray, &numRow);
    double *sol = thomas_solve(triArray, numRow);
    thomas_print(sol, numRow);

    /* Free variables and close files */
    free(triArray);
    free(sol);
    fclose(triSysFile);
    return;
}
/**************************************************************************/

/* Reading in data for tri-diagonal linear system
*/
tridata_t *read_tri_data(FILE *triSysFile, tridata_t* triArray, int *numRow) {

    char ch;
    double val;
    int buffer = INITIALBUFFER, read = 0;

    /* Checking if the data file begins with a number, terminating if so */
    if (fscanf(triSysFile, "%lf", &val) == 1) {
        invalid_data_error("Data formatting error", __LINE__);
    }

    /* Taking the header line from the csv file and throwing it away */
    fscanf(triSysFile, "%c", &ch);
    while (ch != '\n') {
        fscanf(triSysFile, "%c", &ch);
    }

    /* Scan in the rows one at a time and reallocating when needed */
    while (fscanf(triSysFile, "%lf,%lf,%lf,%lf\n", &triArray[read].a,
                                    &triArray[read].b, &triArray[read].c, 
                                    &triArray[read].r) == NUMFIELDSTRI) {
        read++;
        if (read == buffer) {
            buffer *= 2;
            triArray = (tridata_t *)realloc(triArray, sizeof(tridata_t)*buffer);
            assert(triArray != NULL);
        }
    }

    /* Reallocate the mach array to the appropriate size */
    triArray = (tridata_t *)realloc(triArray, sizeof(tridata_t)*read);
    assert(triArray != NULL);

    *numRow = read;
    return triArray;
}
/**************************************************************************/

/* Solving tridiagonal system
*/
double *thomas_solve(tridata_t* triArray, int numRow) {
    int i;
    for (i = 1; i < numRow; i++) {
        triArray[i].b = triArray[i].b-(triArray[i].a*triArray[i-1].c)/triArray[i-1].b;
        triArray[i].r = triArray[i].r-(triArray[i].a*triArray[i-1].r)/triArray[i-1].b;
    }

    double *sol = (double *)malloc(sizeof(double)*numRow);
    assert(sol != NULL);
    sol[numRow-1] = triArray[numRow-1].r/triArray[numRow-1].b;
    for (i = numRow-2; i >= 0; i--) {
        sol[i] = (triArray[i].r-triArray[i].c*sol[i+1])/triArray[i].b;
    }

    return sol;
}
/**************************************************************************/

/* Printing solutions to the tridiagonal system
*/
void thomas_print(double *sol, int numRow) {
    int i;
    FILE *writeFile = safe_fopen("out_linsolve.csv", "w");
    fprintf(writeFile, "x\n");
    for (i = 0; i < numRow; i++) {
        fprintf(writeFile, "%lf\n", sol[i]);
    }
    fclose(writeFile);
}
/**************************************************************************/

/* Function to perform cubic spline interpolation with provided data 
*/
void interp(const char* q4_file, double xo) {
    interpdata_t* interpArray = (interpdata_t *)malloc(sizeof(interpdata_t)*INITIALBUFFER);
    assert(interpArray != NULL);
    FILE *interpFile = safe_fopen(q4_file, "r");
    int numPoints = 0, numIntervals = 0;

    /* Reading in and fitting cubic spline into the data */
    interpArray = read_interp_data(interpFile, interpArray, &numPoints);
    tridata_t* triArray = (tridata_t *)malloc(sizeof(tridata_t)*numPoints);
    assert(triArray != NULL);
    interp_tridiagonal(interpArray, triArray, numPoints);
    double *c = thomas_solve(triArray, numPoints);
    double *b = b_solve(c, interpArray, numPoints);
    double *d = d_solve(c, interpArray, numPoints);
    int *intervals = interp_interval(interpArray, xo, numPoints, &numIntervals);
    interpolate_value(intervals, numIntervals, interpArray, c, b, d, xo);

    /* Free variables and close files */
    free(intervals);
    free(interpArray);
    free(triArray);
    free(c);
    free(b);
    free(d);
    fclose(interpFile);
}
/**************************************************************************/

/* Reading in data for cubic spline interpolation
*/
interpdata_t *read_interp_data(FILE *interpFile, interpdata_t* interpArray, int *numPoints) {

    char ch;
    double val;
    int buffer = INITIALBUFFER, read = 0;

    /* Checking if the data file begins with a number, terminating if so */
    if (fscanf(interpFile, "%lf", &val) == 1) {
        invalid_data_error("Data formatting error", __LINE__);
    }

    /* Taking the header line from the csv file and throwing it away */
    fscanf(interpFile, "%c", &ch);
    while (ch != '\n') {
        fscanf(interpFile, "%c", &ch);
    }

    /* Scan in the rows one at a time and reallocating when needed */
    while (fscanf(interpFile, "%lf,%lf\n", &interpArray[read].x,
                                           &interpArray[read].fx) == NUMFIELDSINTERP) {
        read++;
        if (read == buffer) {
            buffer *= 2;
            interpArray = (interpdata_t *)realloc(interpArray, sizeof(interpdata_t)*buffer);
            assert(interpArray != NULL);
        }
    }

    /* Reallocate the mach array to the appropriate size */
    interpArray = (interpdata_t *)realloc(interpArray, sizeof(interpdata_t)*read);
    assert(interpArray != NULL);

    *numPoints = read;
    return interpArray;
}
/**************************************************************************/

/* Setting up the tridiagonal matrix for the given data
*/
void interp_tridiagonal(interpdata_t* interpArray, tridata_t* triArray, int numPoints) {
    int i;
    for (i = 0; i < numPoints; i++) {
        if (i == 0 || i == numPoints-1) {
            triArray[i].a = 0;
            triArray[i].b = 1;
            triArray[i].c = 0;
            triArray[i].r = 0;
        }
        else {
            triArray[i].a = interpArray[i].x-interpArray[i-1].x;
            triArray[i].c = interpArray[i+1].x-interpArray[i].x;
            triArray[i].b = 2*(triArray[i].a+triArray[i].c);
            triArray[i].r = 3/triArray[i].c*(interpArray[i+1].fx-interpArray[i].fx)+
            3/triArray[i].a*(interpArray[i-1].fx-interpArray[i].fx);
        }
    }
}
/**************************************************************************/

/* Solving the d coefficients for cubic spline
*/
double *d_solve(double* c, interpdata_t* interpArray, int numPoints) {
    int i;
    double *sol = (double *)malloc(sizeof(double)*(numPoints-1));
    assert (sol != NULL);
    for (i = 0; i < numPoints-1; i++) {
        sol[i] = (c[i+1]-c[i])/(3*(interpArray[i+1].x-interpArray[i].x));
    }
    return sol;
}
/**************************************************************************/

/* Solving the b coefficients for cubic spline
*/
double *b_solve(double *c, interpdata_t *interpArray, int numPoints) {
    int i;
    double *sol = (double *)malloc(sizeof(double)*(numPoints-1));
    assert (sol != NULL);
    for (i = 0; i < numPoints-1; i++) {
        sol[i] = (1/(interpArray[i+1].x-interpArray[i].x)*(interpArray[i+1].fx-interpArray[i].fx))-
                 (((interpArray[i+1].x-interpArray[i].x)/3)*(2*c[i]+c[i+1]));
    }
    return sol;
}
/**************************************************************************/

/* Finding the appropriate interval to identify coefficients for interpolation
*/
int *interp_interval(interpdata_t *interpArray, double xo, int numPoints, int *numIntervals) {
    int i, read = 0, buffer = INITIALBUFFER;
    int *interval = (int *)malloc(sizeof(int)*INITIALBUFFER);
    assert(interval != NULL);

    /* Find points and add intervals to dynamic array */
    for (i = 0; i < numPoints-1; i++) {
        if ((interpArray[i].x < xo && interpArray[i+1].x > xo) ||
            (interpArray[i].x > xo && interpArray[i+1].x < xo) ||
            (interpArray[i].x == xo)) {

            interval[read++] = i;
            if (read == buffer) {
                buffer *= 2;
                interval = (int *)realloc(interval, sizeof(int)*buffer);
                assert(interval != NULL);
            }

        }
    }

    /* No interval found */
    if (read == 0) {
        invalid_data_error("Specified interval not in data range", __LINE__); 
    }

    interval = (int *)realloc(interval, sizeof(int)*read);
    assert(interval != NULL);
    *numIntervals = read;
    return interval;
}
/**************************************************************************/

/* Interpolate at given xo using appropriate cubic spline coefficients
*/
void interpolate_value(int *intervals, int numIntervals, interpdata_t *interpArray,
                       double *c, double *b, double *d, double xo) {
    FILE *writeFile = safe_fopen("out_interp.csv", "w");
    int i;
    fprintf(writeFile, "xo,f(xo)\n");
    for (i = 0; i < numIntervals; i++){
        fprintf(writeFile, "%.6lf,%.6lf\n", xo, (interpArray[intervals[i]].fx+b[intervals[i]]*
                                                (xo-interpArray[intervals[i]].x)+
                                                c[intervals[i]]*(xo-interpArray[intervals[i]].x)*
                                                (xo-interpArray[intervals[i]].x)+
                                                d[intervals[i]]*(xo-interpArray[intervals[i]].x)*
                                                (xo-interpArray[intervals[i]].x)*
                                                (xo-interpArray[intervals[i]].x)));
    }
    fclose(writeFile);
    return;
}
/**************************************************************************/

/* Calculate the advection of the fluid using numerical methods
*/
void advection(const char* q5_file) {
    FILE *advecFile = safe_fopen(q5_file, "r");
    advecdata_t advecData = read_advec_data(advecFile);
    double *discreteX = (double *)malloc(sizeof(double)*(advecData.nx+1));
    assert(discreteX != NULL);

    /* Discretise data and calculate advection */
    discretise_intervals(discreteX, advecData);
    calculate_advection(discreteX, advecData);

    /* Free variables and close files */
    free(discreteX);
    fclose(advecFile);
    return;
}
/**************************************************************************/

/* Reading in data for advection data
*/
advecdata_t read_advec_data(FILE *advecFile) {

    char ch;
    double val;
    advecdata_t advecData;

    /* Checking if the data file begins with a number, terminating if so */
    if (fscanf(advecFile, "%lf", &val) == 1) {
        invalid_data_error("Data formatting error", __LINE__);
    }

    /* Taking the header line from the csv file and throwing it away */
    fscanf(advecFile, "%c", &ch);
    while (ch != '\n') {
        fscanf(advecFile, "%c", &ch);
    }

    /* Scan in the advection quantaties */
    if (fscanf(advecFile, "%lf,%d,%lf,%lf", &advecData.u, &advecData.nx, &advecData.CFL,
                                             &advecData.tf) == NUMFIELDSADVEC);
    return advecData;
}
/**************************************************************************/

/* Calculate the time steps and delta x steps for the given data
*/
void discretise_intervals(double *discreteX, advecdata_t advecData) {
    int i;
    for (i = 0; i < (advecData.nx+1); i++){
        if (i == 0){
            discreteX[i] = 0;
        }
        else {
            discreteX[i] = discreteX[i-1] + 1.0/advecData.nx;
        }
    }
    return;
}
/**************************************************************************/

/* Forward Euler and Central Difference schemes to numerically solve the advection equation
*/
void calculate_advection(double *discreteX, advecdata_t advecData) {
    int i, n;
    FILE* writeFile = safe_fopen("out_advection.csv", "w");
    double dx = 1.0/advecData.nx;
    double dt = (advecData.CFL*dx)/advecData.u;
    int Nt = ceil((advecData.u*advecData.tf)/(advecData.CFL*(1.0/(advecData.nx))));
    int Nx = advecData.nx+1;
    double RHS;

    /* Allocating space for phi vectors */
    double *phi = (double *)malloc(Nx*sizeof(double));
    assert(phi != NULL);
    double *prevPhi = (double *)malloc(Nx*sizeof(double));
    assert(prevPhi != NULL);
    double *prevPrevPhi = (double *)malloc(Nx*sizeof(double));
    assert(prevPrevPhi != NULL);


    /* Calculate inital phi conditions */
    for (i = 0; i < Nx; i++) {
        if (discreteX[i] >= 0.125 && discreteX[i] <= 0.375) {
            prevPrevPhi[i] = 0.5*(1-cos(8*M_PI*(discreteX[i]-0.125)));
        }
        else {
            prevPrevPhi[i] = 0;
        }
    }

    /* Calculting system for the first time step with forward Euler scheme */
    for (i = 0; i < Nx; i++) {
        if (i == 0 || i == Nx-1) {
            RHS = -1.0*advecData.u*((prevPrevPhi[1]-prevPrevPhi[Nx-2])/(2*dx));
            prevPhi[i] = RHS*dt + prevPrevPhi[i];
        }
        else {
            RHS = -1.0*advecData.u*((prevPrevPhi[i+1]-prevPrevPhi[i-1])/(2*dx));
            prevPhi[i] = RHS*dt + prevPrevPhi[i];
        }
    }           

    /* Central difference scheme to solve system */
    for (n = 1; n < Nt; n++) {
        for (i = 0; i < Nx; i++) {
            if (i == 0 || i == Nx-1) {
                RHS = -1.0*advecData.u*((prevPhi[1]-prevPhi[Nx-2])/(2*dx));
                phi[i] = RHS*2.0*dt + prevPrevPhi[i];
            }
            else {
                RHS = -1.0*advecData.u*((prevPhi[i+1]-prevPhi[i-1])/(2*dx));
                phi[i] = RHS*2.0*dt + prevPrevPhi[i];
            }
        }

        /* Update phi vectors */
        for (i = 0; i < Nx; i++) {
            prevPrevPhi[i] = prevPhi[i];
            prevPhi[i] = phi[i];
        }
    }

    /* Print out final phi values after Nt time steps */
    fprintf(writeFile, "x,phi\n");
    for (i = 0; i<Nx; i++){
        fprintf(writeFile, "%.6lf,%.6lf\n", discreteX[i], phi[i]);
    }

    /* Free vectors and close files */
    free(phi);
    free(prevPhi);
    free(prevPrevPhi);
    fclose(writeFile);
    return;
}