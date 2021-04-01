/******************************************************************************\
*								 Definitions							 *
\******************************************************************************/
#include <iostream>
#include <iomanip>      
#include <cstdlib>
#include <cmath>
#include <list> 
#define CHAR_LEN 1000
#define EPS1 1.0e-8
#define E  2.7182818284590452353602874713526625
#define INF 1.0e99
#define PI 3.1415926535897932384626433832795029

using namespace std; 

// Data structures
typedef double allele; 										// type allele
typedef struct {
			allele *chromosome;								
			double fitness;									
} individual;												// data structure individual
typedef struct {			
			individual *ind;
			double sum_fitness;
			double mean_fitness;
			double max_fitness;
			int best_individual;		
} population;												// data structure population


// Global variables and parameters
extern int chrom_size;										// size of the chromosome (dimension of the problem) 
extern int gen;												// generation
extern int flag_last_gen;									// flag indicating the last generations
extern int n_instances;										// number of instances 
extern int n_runs_per_instance;								// number of runs for each instance 
extern int popsize;											// size of the population 
extern int save_datagen_flag;								// flag for saving data for generation in the first run
extern int type_crossover;									// crossover type: 1-2PointX; 2-UX; 3-PX, 4-simple ePX, 5-ePX, 6-Multiple 2PointX, 7-Multiple UX, 8-mEX+ePX v1, 9-mEX+ePX v2						
extern long int max_gen;									// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)
extern long int max_rec_comp, tot_rec_comp;					// Statistics: maximum number and total number of recombination components (for PX and rPX)
extern long int randnseed;									// seed for Gaussian rand generator
extern long int sum_sucRate, sum_impRate, sum_cross;		// Statistics: sum of sucessful recombinations, improvements, crossovers
extern double F_DE;											// differential weight: DE parameter
extern population popc;										// current population 
// Vectors
extern int *vsort_aux;										// auxiliar sorted vector of integers (used in different methods)
extern int *file_gen;										// data to be stored: number of generations								
extern double *file_best_fitness, *time_run;				// data to be stored: best fitness, runtime
extern double *file_best_fitness_gen;						// data to be stored: best fitness over the generations for run 0
extern double *file_impRate, *file_sucRate;					// data to be stored: improvement and successful crossover rates 
// Matrices
extern double **File_best_ind;								// data to be stored: best individual
extern double **File_recComp;								// data to be stored: recombination partition information
// CEC17 Test Function Suite for Single Objective Optimization
extern int ini_flag,n_flag,func_flag,*SS;
extern int func_num;
extern double *OShift,*M,*y_cec2017,*z_cec2017,*x_bound;

// Function declaration
//  aux_functions.cpp
void desaloc_matrixd(double **Matrix , int lines);
void desaloc_matrixi(int **Matrix , int lines);
void rand_perm_size(int *inp, int *out, int size_inp, int size_out);
void rotMatrix1(double **A, int l);
int isEqual(double x, double y);
int isVectorEqual(int *v1, int *v2, int size_v);
int random_int(int L_range, int H_range);
int *aloc_vectori(int lines);
int **aloc_matrixi(int lines , int collums);
double gasdev(long int *idum);
double random_dou(void);
double *aloc_vectord(int lines);
double **aloc_matrixd(int lines , int collums);
individual *aloc_vectorind(int lines);
// file_man.cpp
void file_output(int N_Mk, int model_Mk, int total_runs);
// statistics.cpp
void statistics(population *pop, int n_run);
// transformation.cpp
void Point2X(allele *parent1, allele *parent2, allele *offspring);
void UX(allele *parent1, allele *parent2, allele *offspring);
// cec17_test_func.cpp
void cec17_test_func_init(int nx, int mx);
void cec17_test_func(double *x, double *f, int nx, int mx);

