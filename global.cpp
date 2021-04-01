#include "defs.h"

// Global variables and parameters
int chrom_size;										// size of the chromosome (dimension of the problem) 
int flag_last_gen;									// flag indicating the last generations
int gen;											// generation
int n_instances = 1;								// number of instances 
int n_runs_per_instance = 50;						// number of runs for each instance 
int popsize=100;									// size of the population 
int save_datagen_flag=0;							// flag for saving data for generation in the first run
int type_crossover;									// crossover type: 1-2PointX; 2-UX; 3-PX, 4-simple ePX, 5-ePX, 6-Multiple 2PointX, 7-Multiple UX, 8-mEX+ePX v1, 9-mEX+ePX v2	
long int max_gen;									// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)
long int max_rec_comp, tot_rec_comp;				// Statistics: maximum number and total number of recombination components (for PX and rPX)
long int randnseed;									// seed for Gaussian rand generator
long int sum_sucRate, sum_impRate, sum_cross;		// Statistics: sum of sucessful recombinations, improvements, crossovers
double F_DE=0.8;									// differential weight: DE parameter
population popc;									// current population 
// Vectors
int *vsort_aux;										// auxiliar sorted vector of integers (used in different methods)
int *file_gen;										// data to be stored: number of generations								
double *file_best_fitness, *time_run;				// data to be stored: best fitness, runtime
double *file_best_fitness_gen;						// data to be stored: best fitness over the generations for run 0
double *file_impRate, *file_sucRate;				// data to be stored: improvement and successful crossover rates 
// Matrices
double **File_best_ind;								// data to be stored: best individual
double **File_recComp;								// data to be stored: recombination partition information
// CEC17 Test Function Suite for Single Objective Optimization
int ini_flag=0,n_flag,func_flag,*SS;
int func_num;
double *OShift,*M,*y_cec2017,*z_cec2017,*x_bound;
