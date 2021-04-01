*** Standard Differential Evolution  with epsilon-Partition Crossover (ePX) ***

Description: This is the source code for the DE with ePX for bound constrained optimization functions of the CEC17 Suite for Single Objective Optimization. 

Reference:  Tinos, R.; Whitley, D.; Chicano, F. & Ochoa, G. (2021), "Partition Crossover for Continuous Optimization: ePX", Proc. of GECCO'2021.	

Contact: Renato Tinos <rtinos at ffclrp.usp.br>


Running the code: 

	./de_rpx_bc <crossover> <function> <N> (optional: <model_Mk>)"<<endl

<crossover>: type of crossover. 1-2PointX; 2-UX; 3-PX, 4-simple ePX (without removing common vertices), 5-ePX, 6-Multiple 2PointX, 7-Multiple UX.	

<function>: CEC17 Test Function number. 1-Bent Cigar; 4-Rosenbrock's Function; 5-Rastrigin's Function; 10-Schwefel's Function.		    

<N>: number of dimensions. 10; 20; 30; 50; 100.

<K>: controls the epistasis degree of the NK landscapes instance. K is a non-negative integer smaller than or equal to N.

<model_NK>: epistasis types. 0-no rotation (default); 1-rotation: for the CEC17 Test Function Suite; 2-rotation: adjacent with sobreposition. Rotations 1 and 2 are not allowed in function 4.


Example for running the code for: Rosenbrock's Function (f4) with N=30, no rotation, and crossover ePX:

make

./de_ePX_bc 5 4 30 0


Observation 1: In order to run the DE for the CEC17 Suite for Single Objective Optimization, it is necessary to copy the directory "\input_data" from codes.rar at https://github.com/P-N-Suganthan/CEC2017-BoundContrained


Observation 2: Class Mk, given in cMk.h, implements Mk landscapes (that is a generalization of k-bounded functions). Some functions in Mk.h:

- double Mk::ePX(double *parent1, double *parent2, double *offspring, double alpha);  epsilon-Partition Crossover. Inputs: two parents, pointer to offspring, and alpha=1-epsilon. Outputs: offspring and its fitness.
	
- double Mk::PX(double *parent1, double *parent2, double *offspring);  Partition Crossover. Inputs: two parents, pointer to offspring. Outputs: offspring and its fitness.

- double Mk::compFitness (double *x); Fitness computation for Mk landscapes representing k-bounded single objective functions. Here, it calls cec17_compFitness(.), of cec17_subfunc.h, that is adapted from the CEC17 Test Function Suite for Single Objective Optimization. Input: solution x; Output: fitness of x. 
		

Observation 3: file global.cpp contains the parameters of the DE (examples: population size and parameter F) and problem (example: number of runs). Parameter alpha (i.e, 1-epsilon) is given in function crossover() of de.cpp.


Observation 4: de_ePX_bc generates 6 main files
 
- bfi_f%d_N%d_m%d_c%d.dat (example for ./de_ePX_bc 5 4 30 1: bfi_f4_N30_m1_c5.dat): best fitness found in each run
	
- bind_f%d_N%d_m%d_c%d.dat: best individuals found in each run

- time_f%d_N%d_m%d_c%d.dat: time for each run

- sucRate_f%d_N%d_m%d_c%d.dat: successful recombination rate for each run

- impRate_f%d_N%d_m%d_c%d.dat: improvement rate for each run

- recComp_f%d_N%d_m%d_c%d.dat: saves the number of recombination components for PX, ePX and sePX