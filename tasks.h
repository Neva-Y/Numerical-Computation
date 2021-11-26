/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : 1001969
 *   Name        : Goh Xing Yang
 *
 ***************************************************************************/

/* Header guards
*/
#ifndef TASKS_H
#define TASKS_H

/************************************************************************/

/* Hash defined constants
*/
#define SHOCKFILE 1
#define LINSOLVEFILE 2
#define INTERPFILE 3
#define INTERPPOINT 4
#define ADVECTIONFILE 5
#define M_PI 3.14159265358979323846
#define PI_DEGREES 180
#define MAXITER 10000
#define EPS 1e-10
#define NUMFIELDSSHOCK 5
#define NUMFIELDSTRI 4
#define NUMFIELDSINTERP 2
#define NUMFIELDSADVEC 4
#define INITIALBUFFER 5
#define DETATCHED -1
#define THETAMAX 90
#define THETAMIN 0

/************************************************************************/

/* Struct definition for input data fields
*/
typedef struct {
    double M;
    double theta;
    double beta_l;
    double beta_u;
    double gamma;
} machdata_t;

typedef struct {
	double a;
	double b;
	double c;
	double r;
} tridata_t;

typedef struct {
	double x;
	double fx;
} interpdata_t;

typedef struct {
	double u;
	int nx;
	double CFL;
	double tf;
} advecdata_t;

/************************************************************************/

/* Function prototypes to open files and print errors 
*/
FILE *safe_fopen(const char* path, const char* mode);
void invalid_data_error(char *msg, int line);

/* Shockwave functions 
*/
double f_beta(double M, double theta, double gamma, double beta);
double f_beta_prime(double M, double theta, double gamma, double beta);
void read_shockwave_a(FILE *shockFile, machdata_t* dataLine);
double *read_mach_array(FILE *shockFile, double *machArray, int *numMach);
double newtonraph(double guess, double theta, double M, double gamma);
void print_shockwave(machdata_t* dataLine, double *machArray, int numMach);
void shockwave(const char* q1_file);

/* Tri-diagonal linear solving functions 
*/
tridata_t *read_tri_data(FILE *triSysFile, tridata_t* triArray, int *numRow);
double *thomas_solve(tridata_t* triArray, int numRow);
void thomas_print(double *sol, int numRow);
void linsolve(const char* q3_file);

/* Cubic spline interpolation functions 
*/
interpdata_t *read_interp_data(FILE *interpFile, interpdata_t* interpArray, int *numPoints);
void interp_tridiagonal(interpdata_t* interpArray, tridata_t* triArray, int numPoints);
double *b_solve(double* c, interpdata_t* interpArray, int numPoints);
double *d_solve(double* c, interpdata_t* interpArray, int numPoints);
int *interp_interval(interpdata_t* interpArray, double xo, int numPoints, int *numIntervals);
void interpolate_value(int *intervals, int numIntervals, interpdata_t *interpArray,
                       double *c, double *b, double *d, double xo);
void interp(const char* q4_file, double xo);

/* Advection equation solving functions 
*/
advecdata_t read_advec_data(FILE *advecFile);
void calculate_advection(double *discreteX, advecdata_t advecData);
void discretise_intervals(double *discreteX, advecdata_t advecData);
void advection(const char* q5_file);

#endif