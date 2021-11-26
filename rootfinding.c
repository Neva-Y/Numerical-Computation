#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>

/**************************************************************************/

#define SHOCKFILE 1
#define LINSOLVEFILE 2
#define INTERPFILE 3
#define ADVECTIONFILE 4
#define M_PI 3.14159265358979323846
#define MAXITER 10000
#define EPS 1e-6
#define NUMFIELDS 5
#define MACHBUFFER 5
#define DETATCHED -1
#define THETAMAX 90
#define THETAMIN 0

/**************************************************************************/

/* Defining the data fields for the first section of the input file 
*/
typedef struct {
    double M;
    double theta;
    double beta_l;
    double beta_u;
    double gamma;
} data_t;

/**************************************************************************/

FILE* safe_fopen(const char* path, const char* mode);
double f_beta(double M, double theta, double gamma, double beta);
double f_beta_prime(double M, double theta, double gamma, double beta);
void invalid_data_error(char *msg, int line);
void read_shockwave_a(FILE *shockFile, data_t* dataLine);
void switchArrayDoubles(double *machArray, int i, int j);
double *read_mach_array(FILE *shockFile, double *machArray, int *numMach);
double newtonraph(double guess, double theta, double M, double gamma);
void print_shockwave(data_t* dataLine, double *machArray, int numMach);

/**************************************************************************/

/* Main binding function
*/
int main(int argc, char *argv[]) {
	data_t* dataLine = (data_t *)malloc(sizeof(data_t));
	assert(dataLine != NULL);
    FILE *shockFile = safe_fopen(argv[SHOCKFILE], "r");
	int numMach = 0;
	double *machArray = (double *)malloc(sizeof(double)*MACHBUFFER);
	assert(machArray != NULL);

	read_shockwave_a(shockFile, dataLine);
	machArray = read_mach_array(shockFile, machArray,&numMach);

    printf("%lf,%lf,%lf,%lf,%lf", dataLine->M, dataLine->theta,
                                    dataLine->beta_u, dataLine->beta_l,
                                    dataLine->gamma);

    printf("\n%d",numMach);
    for (int i=0; i<numMach; i++) {
        printf("\n%lf\n", machArray[i]);
    }
    print_shockwave(dataLine, machArray, numMach);
    free(machArray);
    free(dataLine);
	return 0;
}
/**************************************************************************/

/* Function to safely open a file
*/
FILE* safe_fopen(const char* path, const char* mode) {
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
the first section of the data
*/
void read_shockwave_a(FILE *shockFile, data_t* dataLine) {

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
    if (fscanf(shockFile, "%lf,%lf,%lf,%lf,%lf\n", &dataLine->M, &dataLine->theta,
                                    &dataLine->beta_u, &dataLine->beta_l,
   	                                &dataLine->gamma)!=5) {
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

/* Dynamically read the Mach numbers for section b
*/
double *read_mach_array(FILE *shockFile, double *machArray, int *numMach) {
	int i, j, buffer = MACHBUFFER, read = 0;

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

	/* Sorting the mach numbers in ascending order */
	for(i=1; i<read; i++){
        j=i;
        while(j>0 && machArray[j] < machArray[j-1]){
            switchArrayDoubles(machArray, j, j-1);
            j--;
        }
    }
    *numMach = read;
	return machArray;
}
/**************************************************************************/

/* Function to swap two doubles with pointers
*/
void switchArrayDoubles(double *machArray, int i, int j) {
    double tmp=machArray[i];
    machArray[i]=machArray[j];
    machArray[j]=tmp;
    return;
}
/**************************************************************************/

double f_beta(double M, double theta, double gamma, double beta) {
    return (2*(1/tan(beta))*((M*M*(sin(beta))*(sin(beta))-1)/(M*M*(gamma + cos(2*beta))+2))-tan(theta));
}
/**************************************************************************/

double f_beta_prime(double M, double theta, double gamma, double beta) {
    return (4*M*M*cos(beta)*sin(beta))/(tan(beta)*((gamma + cos(2*beta))*M*M + 2)) - (2*(tan(beta)*tan(beta) + 1)*(M*M*sin(beta)*sin(beta) - 1))/(tan(beta)*tan(beta)*((gamma + cos(2*beta))*M*M + 2)) + (4*M*M*sin(2*beta)*(M*M*sin(beta)*sin(beta) - 1))/(tan(beta)*(M*M*(gamma + cos(2*beta)) + 2)*(M*M*(gamma + cos(2*beta)) + 2));

}
/**************************************************************************/

double newtonraph(double guess, double theta, double M, double gamma) {
    double fx1, fpx1;
    int iterations=0;
    double beta = guess*M_PI/180;
    fx1 = f_beta(M, theta*M_PI/180, gamma, beta); 
    fpx1= f_beta_prime(M, theta*M_PI/180, gamma, beta);
    
    while (fabs(fx1) > EPS && iterations < MAXITER) {
        iterations++;
        beta = beta-(fx1/fpx1);
        fx1 = f_beta(M, theta*M_PI/180, gamma, beta);
        fpx1= f_beta_prime(M, theta*M_PI/180, gamma, beta);
    }
    //printf("Iternations: %d \n",iterations);
    if (iterations == MAXITER && fabs(fx1) > EPS) {
        return DETATCHED;
    }
    return beta*180/M_PI;
}
/**************************************************************************/

void print_shockwave(data_t* dataLine, double *machArray, int numMach) {
    FILE *writeFile = safe_fopen("out_shock.csv", "w");
    int i;
    double bl, bu, theta;
    fprintf(writeFile, "M,theta,beta_lower,beta_upper\n");
    for (i = 0; i < numMach; i++) {

        /* Initial guesses for theta = 0 */
        theta = 0;
        bl = asin(1/machArray[i])*180/M_PI;
        bu = THETAMAX;

        /* Calculate angles for increments of theta */
        bl = newtonraph(bl, theta, machArray[i], dataLine->gamma);
        bu = newtonraph(bu, theta, machArray[i], dataLine->gamma);

        while (bl!=DETATCHED && bu!=DETATCHED) {
            fprintf(writeFile, "%0.4lf,%0.4lf,%0.4lf,%0.4lf\n", machArray[i], theta, bl, bu);
            theta++;

            bl = newtonraph(bl, theta, machArray[i], dataLine->gamma);
            bu = newtonraph(bu, theta, machArray[i], dataLine->gamma);

            if (bl<THETAMIN || bu>THETAMAX || bu<bl) {
                break;
            }

        }

        if (i < numMach-1) {
            fprintf(writeFile, ".\n.\n");   
        }
    }
    fclose(writeFile);
    return;
}