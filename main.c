/***************************************************************************
 *
 *   File        : main.c
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

int main(int argc, char *argv[]) {
	
	/* Parsing in required information from command line */
	char* q1_file = argv[SHOCKFILE];
	char* q3_file = argv[LINSOLVEFILE];
	char* q4_file = argv[INTERPFILE];
	double xo = atof(argv[INTERPPOINT]);
	char* q5_file = argv[ADVECTIONFILE];

	/* Initialising structs defined in sys/time.h to record elapsed time */
    struct timeval start;
    struct timeval stop;  
    
	/* Question 1 */
	gettimeofday(&start, NULL);
	shockwave(q1_file);
	gettimeofday(&stop, NULL);
    double elapsedMs = (stop.tv_sec - start.tv_sec) * 1000.0;
    elapsedMs += (stop.tv_usec - start.tv_usec) / 1000.0;
    printf("Question 1 Shockwave Equation:  %.2f milliseconds\n", elapsedMs);

	/* Question 3 */
	gettimeofday(&start, NULL);
	linsolve(q3_file);
	gettimeofday(&stop, NULL);
	elapsedMs = (stop.tv_sec - start.tv_sec) * 1000.0;
    elapsedMs += (stop.tv_usec - start.tv_usec) / 1000.0;
    printf("Question 3 Tri-diagonal Linear System:  %.2f milliseconds\n", elapsedMs);

	/* Question 4 */
	gettimeofday(&start, NULL);
	interp(q4_file, xo);
	gettimeofday(&stop, NULL);
    elapsedMs = (stop.tv_sec - start.tv_sec) * 1000.0;
    elapsedMs += (stop.tv_usec - start.tv_usec) / 1000.0;
    printf("Question 4 Cubic Spline Interpolation:  %.2f milliseconds\n", elapsedMs);

	/* Question 5 */
	gettimeofday(&start, NULL);
	advection(q5_file);
	gettimeofday(&stop, NULL);
    elapsedMs = (stop.tv_sec - start.tv_sec) * 1000.0;
    elapsedMs += (stop.tv_usec - start.tv_usec) / 1000.0;
    printf("Question 5 Advection Equation:  %.2f milliseconds\n", elapsedMs);
    
	return (EXIT_SUCCESS);
}
