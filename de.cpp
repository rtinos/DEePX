/*******************************************************************************\
*  	Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 		  	   	*
*																				*
* 	Standard Differential Evolution  with epsilon-Partition Crossover (ePX)		*
* 	Test problem: bound constrained optimization functions of CEC17 Suite for   *
*					Single Objective Optimization								*
*						 														*
* 	Copyright (C) 2021  Renato Tinos <rtinos@ffclrp.usp.br>						*
* 						 														*
* 	Reference: Tinos, R.; Whitley, D.; Chicano, F. & Ochoa, G. (2021), 		 	*                    
* 				"Partition Crossover for Continuous Optimization: ePX",		 	*
*				Proc. of GECCO'2021.											*
*																				*
* 	de_epx_bc is free software: you can redistribute it and/or modify it 		*
* 		under the terms of the GNU General Public License as published by the	*
* 		Free Software Foundation, either version 3 of the License, or			*
* 		(at your option) any later version.						 				*
* 						 														*
* 	de_epx_bc is distributed in the hope that it will be useful, but			*
* 		WITHOUT ANY WARRANTY; without even the implied warranty of				*
* 		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.					*
* 		See the GNU General Public License for more details.					*
* 																				*
* 	You should have received a copy of the GNU General Public License along		*
* 		with this program.  If not, see <http://www.gnu.org/licenses/>.			*		
\*******************************************************************************/

#include <time.h>
#include "defs.h"
#include "cMk.h"


/******************************************************************************\
*				  	Print population					 			 			*
\******************************************************************************/
void print_data(population *pop, int n_run){

	cout <<"Generation:"<< gen << ", run: "<<n_run<<endl;
	cout <<"Best individual:"<< pop->best_individual << endl;
	cout <<"Fitness of the best individual:"<< pop->max_fitness << endl;
	cout <<"Mean fitness: "<< pop->mean_fitness << endl;	
	/*for (int i=0;i<popsize ;i++) {	
		cout <<"("<< pop->ind[i].fitness<<") " ;
		for (int gene=0;gene<chrom_size ;gene++) 
			cout << pop->ind[i].chromosome[gene]<<" ";
		cout << endl;
	}*/
	
}


/******************************************************************************\
*								Fitness Computation    				   		   *
\******************************************************************************/

double compFitness( allele *ind ){
	double f[1], fitness;
		
	cec17_test_func(ind, f, chrom_size ,1);
	fitness=-f[0];			// transforming to a maximization problem
	
	return fitness;
}

/***********************************************************************************************************\
*								 Crossover	 													     		*
* Crossover type: 1-2PointX; 2-UX; 3-PX, 4-simple ePX, 5-ePX, 6-Multiple 2PointX, 7-Multiple UX			    *
*                 8-mEX+ePX v1, 9-mEX+ePX v2								   								*
\***********************************************************************************************************/
double crossover(Mk *Mk_instance, allele *parent1, allele *parent2, allele *offspring)
{
	int n_rep=10; // number of repetitions for types 6 and 7
	int type_crossover_aux;
	double fit_offspring, fit_offspring_aux, *offspring_aux, alpha=0.1; 

	if (type_crossover==8){
		// mEX+ePX v1
		if (flag_last_gen==1)
			type_crossover_aux=5;
		else
			type_crossover_aux=6;
	}
	else if (type_crossover==9){
	    // mEX+ePX v2	
		if (random_dou()<=0.05)
			type_crossover_aux=5;
		else
			type_crossover_aux=6;				
	}
	else {
		type_crossover_aux=type_crossover;
	}

	if (type_crossover_aux==1){
		// 2-point Crossover
		Point2X(parent1,parent2,offspring);
		fit_offspring=compFitness(offspring);	
	}
	else if (type_crossover_aux==2){
		// Uniform Crossover
		UX(parent1,parent2,offspring);
		fit_offspring=compFitness(offspring);	
	}	
	else if (type_crossover_aux==3){
		// Partition Crossover
		fit_offspring=Mk_instance->PX(parent1,parent2,offspring);	
	}
	else if (type_crossover_aux==4){
		// simple epsilon-Partition Crossover (i.e., without removing common vertices)
		fit_offspring=Mk_instance->sePX(parent1,parent2,offspring,alpha);	
	}
	else if (type_crossover_aux==5){
		// epsilon-Partition Crossover
		fit_offspring=Mk_instance->ePX(parent1,parent2,offspring,alpha);	
	}
	else if (type_crossover_aux==6){
		// Multiple 2-point Crossover
		offspring_aux=aloc_vectord(chrom_size);
		for (int i=0;i<n_rep;i++){
			if (i==0){
				Point2X(parent1,parent2,offspring);
				fit_offspring=compFitness(offspring);
			}
			else{
				Point2X(parent1,parent2,offspring_aux);
				fit_offspring_aux=compFitness(offspring_aux);
				if (i==0 || (fit_offspring_aux>fit_offspring)){
					fit_offspring=fit_offspring_aux;
					for(int gene=0;gene<chrom_size;gene++)
						offspring[gene]=offspring_aux[gene];						
				}
			}
		}
		delete [] offspring_aux;	
	}
	else if (type_crossover_aux==7){
		// Multiple Uniform Crossover
		offspring_aux=aloc_vectord(chrom_size);
		for (int i=0;i<n_rep;i++){
			if (i==0){
				UX(parent1,parent2,offspring);
				fit_offspring=compFitness(offspring);
			}
			else{
				UX(parent1,parent2,offspring_aux);
				fit_offspring_aux=compFitness(offspring_aux);
				if (i==0 || fit_offspring_aux>fit_offspring ){
					fit_offspring=fit_offspring_aux;
					for(int gene=0;gene<chrom_size;gene++)
						offspring[gene]=offspring_aux[gene];						
				}
			}
		}
		delete [] offspring_aux;	
	}		
	
	if ( type_crossover==3 || type_crossover==4 || type_crossover==5 ){
		if (Mk_instance->n_rec_comp > max_rec_comp)
			max_rec_comp=Mk_instance->n_rec_comp;
		tot_rec_comp+=Mk_instance->n_rec_comp;
	}
	
	return fit_offspring;
	
}


/******************************************************************************\
*								 Generation of the DE						   *
\******************************************************************************/
void generation(Mk *Mk_instance, int n_run){
	int j=0 , a, b, c, popsize_1;
	individual y, xnew;
	
	y.chromosome = aloc_vectord(chrom_size);
	xnew.chromosome = aloc_vectord(chrom_size);
	
	popsize_1=popsize-1;
	do {	
		
		// Selection of the Parents
		a=j;
		while(a==j)		
			a=random_int(0, popsize_1);
		b=j;
		while(b==j || b==a)		
			b=random_int(0, popsize_1);
		c=j;
		while(c==j || c==a || c==b)		
			c=random_int(0, popsize_1);
		
		// Transformation
		for (int gene=0;gene<chrom_size;gene++)
			y.chromosome[gene] = popc.ind[a].chromosome[gene] + F_DE*(popc.ind[b].chromosome[gene]-popc.ind[c].chromosome[gene]);		
		y.fitness=compFitness(y.chromosome);
		xnew.fitness=crossover(Mk_instance, y.chromosome, popc.ind[j].chromosome, xnew.chromosome);	
													
		// Statistics: crossover rates
		sum_cross++;	// total number of recombinations
		if ( (xnew.fitness-y.fitness)>EPS1 && (xnew.fitness-popc.ind[j].fitness)>EPS1 )
			sum_sucRate++;	// Record successful recombination rate
		if ( (xnew.fitness-file_best_fitness[n_run])>EPS1 )	
			sum_impRate++; // Record improvement rate
		
		// Selection of Offspring
		if (xnew.fitness>popc.ind[j].fitness){
			for (int gene=0;gene<chrom_size;gene++)
				popc.ind[j].chromosome[gene] = xnew.chromosome[gene];
			popc.ind[j].fitness = xnew.fitness;		
		}
		
		// check if improved best individual
		if ( popc.ind[j].fitness>file_best_fitness[n_run]){		
			file_best_fitness[n_run] = popc.ind[j].fitness;
			for(int i=0;i<chrom_size;i++)
				File_best_ind[n_run][i]=popc.ind[j].chromosome[i];
		}
				
		j = j + 1;	
	} while ( j < popsize);
	
	delete [] y.chromosome;
	delete [] xnew.chromosome;

}	


/******************************************************************************\
*				  	Random individual					 				 	   *
\******************************************************************************/
void randomInd(int num_ind){ 
	double lim_sup, lim_inf;
	
	lim_sup=x_bound[0];
	lim_inf=-x_bound[0];
	
	for (int gene=0;gene<chrom_size;gene++)
     	popc.ind[num_ind].chromosome[gene] = (lim_sup-lim_inf)*random_dou()+lim_inf;		// x_bound is defined in void cec17_test_func_init(int nx, int mx)
    popc.ind[num_ind].fitness = compFitness( popc.ind[num_ind].chromosome );											
	
}


/******************************************************************************\
*				  	Initiate Population 					 				 *
\******************************************************************************/
void initiatePop(int n_run){
				
	// Dynamic allocation: current population
	popc.ind = aloc_vectorind(popsize);

	for (int num_ind=0;num_ind<popsize;num_ind++){
		// Dynamic allocation: chromosomes	
		popc.ind[num_ind].chromosome = aloc_vectord(chrom_size);

		// Random Initialization
		randomInd(num_ind);	 	      				
	}
	file_best_fitness[n_run]=popc.ind[0].fitness;
	for(int i=0;i<chrom_size;i++)
     	File_best_ind[n_run][i]=popc.ind[0].chromosome[i];
	statistics(&popc, n_run);
	//print_data(&popc, n_run);
	
}


/******************************************************************************\
*				  	Run of the DE 			 								   *
\******************************************************************************/
void de(Mk *Mk_instance, int n_run ){
	
	double time_aux; 
	clock_t time_start;
	
	// Statistics of the number of partitions
	tot_rec_comp=0;
	max_rec_comp=0;
	// Statistics: crossover
	sum_sucRate=0;
	sum_impRate=0;
	sum_cross=0;
		
	// Initialization
	time_start=clock();	
	gen=0;
	flag_last_gen=0;
	initiatePop(n_run);								// initiating population
	
	do {
		gen++; 										// generation index
			
		generation(Mk_instance,n_run);
					
		statistics(&popc,n_run);	
		//print_data(&popc,n_run);
					
		time_aux = ( (double) ( clock() - time_start ) ) / ( (double) CLOCKS_PER_SEC);
		if (time_aux >= 0.95*chrom_size) 
			flag_last_gen=1;
	} while ( time_aux < chrom_size );
	//} while ( gen < max_gen );

	
	//time_aux = ( (double) ( clock() - time_start ) ) / ( (double) CLOCKS_PER_SEC);
	// Data to be saved
	time_run[n_run]=time_aux;
	file_gen[n_run]=gen;
	if (gen<max_gen-1){
		// Save data for crossover statistics (if not saved function statistics() yet )
		if (sum_cross>0){
			file_sucRate[n_run]=((double) sum_sucRate)/sum_cross;
			file_impRate[n_run]=((double) sum_impRate)/sum_cross;
		}
		else{
			file_sucRate[n_run]=0.0;
			file_impRate[n_run]=0.0;		
		}
		if (type_crossover==3 || type_crossover==4 || type_crossover==5){
			if (sum_cross>0)
				File_recComp[n_run][0]=((double) tot_rec_comp)/sum_cross;
			else
				File_recComp[n_run][0]=0.0;
			File_recComp[n_run][1]=(double) max_rec_comp;
		}
	}
	
	// Deleting population
	for (int num_ind=0;num_ind<popsize;num_ind++)
		delete [] popc.ind[num_ind].chromosome;	
	delete [] popc.ind;
	
}


/******************************************************************************\
*				  	Mk_model  												   *
*  model_neig 0: no rotation (standard)						   				   *
*  model_neig 1: rotation - CEC17						   				   	   *
*  model_neig 2: rotation - adjacent with sobreposition						   *
\******************************************************************************/
void Mk_model(int N, int model_neig, list<int> *PHI) 
{
	int j, l; 
	double **M_aux;
	
	// Modify rotation matrix M 
	if (model_neig==0){
		// no linear transformation (standard)
		for (l=0;l<N;l++)
			for (j=0;j<N;j++)
				if (l==j)
					M[l*N+l]=1.0;
				else
					M[l*N+j]=0.0;

	}
	else if (model_neig==2){
		// rotation: adjacent 
		M_aux=aloc_matrixd(N,N);		
		rotMatrix1(M_aux,N);
		for (l=0;l<N;l++)
			for (j=0;j<N;j++)
			   M[l*N+j]=M_aux[l][j];
		desaloc_matrixd(M_aux,N);
	}	
	// creates vector of lists PHI (PHI[l] indicates the variables that influences subfunction fl) 
	if (func_num==4){
		for (l=0;l<N-1;l++){
			PHI[l].push_back(l); 						// variable l influences f_l
			PHI[l].push_back(l+1); 						// variable l+1 influences f_l
		}
	}
	else{
		for (l=0;l<N;l++)
			for (j=0;j<N;j++)
				if ( isEqual(M[l*N+j],0.0)==0 )
					PHI[l].push_back(j); 			// variable j influences f_l
	}
		
	/*cout<<"M:"<<endl;
	cout<<fixed;
	for (l=0;l<N;l++) {
		for (j=0;j<N;j++) 
			cout<<M[l*N+j]<<"\t";
		cout<<endl;
	}	
	cout<<"Subfunction: variables "<<endl;
	for (l=0;l<N;l++){
		cout<<" f_"<<l<<": ";
		for(list<int>::iterator jj = PHI[l].begin(); jj != PHI[l].end(); jj++) 
			cout<<"  x_"<<*jj;
		cout<<endl;
	}*/
		
}


/******************************************************************************\
*				  	Main													   *
\******************************************************************************/
int main(int argc , char *argv[])
{
	int total_runs, n_run=0;
	int N_Mk, model_Mk;					// variables for the Mk model
	list<int> *PHI_Mk;					// 	vector of lists PHI (PHI[l] indicates the variables that influences subfunction fl)
	
	
	// Arguments
	if( argc!=4 && argc!=5) {
		cout<<"Incorret number of arguments!"<<endl;
		cout<<"Call: de_rpx_bc <crossover> <function> <N> (optional: <model_Mk>)"<<endl;
		exit(1);
	}
	else{
		type_crossover=atoi(argv[1]);
		if ( type_crossover<1 || type_crossover>9 ) {
			cout<<"Call: de_epx_bc <crossover> <function> <N> (optional: <model_Mk>)"<<endl;
			cout<<"Incorrect argument!"<<endl;
			cout<<"<crossover> types:  1-2PointX; 2-UX; 3-PX, 4-simple ePX, 5-ePX, 6-Multiple 2PointX, 7-Multiple UX, 8-mEX+ePX v1, 9-mEX+ePX v2"<<endl;
			exit(1);
		}
		func_num=atoi(argv[2]);
		if ( func_num!=1 && func_num!=4 && func_num!=5 && func_num!=10  ) {
			cout<<"Call: de_epx_bc <crossover> <function> <N> (optional: <model_Mk>)"<<endl;
			cout<<"Incorrect argument!"<<endl;
			cout<<"<function> CEC17 Test Function number: 1-Bent Cigar;";
			cout<<" 4-Rosenbrock's Function; ";
			cout<<" 5-Rastrigin's Function; ";
		    cout<<" 10-Schwefel's Function; ";		    
			cout<<endl;
			exit(1);
		}	
		N_Mk=atoi(argv[3]);
		if ( N_Mk!=10 && N_Mk!=20 && N_Mk!=30 && N_Mk!=50 && N_Mk!=100 ) {
			cout<<"Call: de_epx_bc <crossover> <function> <N> (optional: <model_NK>)"<<endl;
			cout<<"Incorrect argument!"<<endl;
			cout<<"<N> dimensions: 10; 20; 30; 50; 100"<<endl;
			exit(1);
		}
		if( argc==5) {
			model_Mk=atoi(argv[4]);	
			if (  func_num==4 && model_Mk!=0  ){
				cout<<"Call: de_epx_bc <crossover> <function> <N> (optional: <model_NK>)"<<endl;
				cout<<"Incorrect argument!"<<endl;
				cout<<"function number 4 does not accept parameter <model_NK> different from 0"<<endl;
				exit(1);				
			}
			if  (model_Mk <0 || model_Mk>2 ) {
				cout<<"Call: de_epx_bc <crossover> <function> <N> (optional: <model_NK>)"<<endl;
				cout<<"Incorrect argument!"<<endl;
				cout<<"<model_Mk> epistasis types: 0-no rotation (default); 1-rotation: for CEC17 Test Function Suite; 2-rotation: adjacent with sobreposition"<<endl;
				exit(1);
			}										
		}
		else{
			model_Mk=0;
		}
	}	
		
	// Parameters
	chrom_size=N_Mk;													// size of the chromosome
	total_runs=n_instances*n_runs_per_instance;							// number of runs
	max_gen=(long int)(chrom_size*(10000.0/popsize));					// maximum number of generations (when stop criterion is time, it is used for saving some statistics data)	
																		// 		In the competition CEC17, the number of evaluations is 10000*N
	// Allocation of vectors and matrices
	file_best_fitness_gen=aloc_vectord(max_gen);
	file_best_fitness=aloc_vectord(total_runs);
	file_gen=aloc_vectori(total_runs);
	time_run=aloc_vectord(total_runs);
	file_sucRate=aloc_vectord(total_runs);
	file_impRate=aloc_vectord(total_runs);
	File_best_ind=aloc_matrixd(total_runs,chrom_size);
	File_recComp=aloc_matrixd(total_runs,2);
	vsort_aux=aloc_vectori(chrom_size);									// Auxiliar sorted vector of integers (used in different methods)
	for (int i=0;i<chrom_size;i++)
		vsort_aux[i]=i;

	// Initialization of CEC17 Test Functions
	cec17_test_func_init(chrom_size, 1);

	cout << "\n ***** Standard Differential Evolution ****" << endl;
	cout << "Function "<<func_num<<", N="<<N_Mk<<endl;
	cout << " model_Mk="<<model_Mk<< endl;
	cout << "type_crossover="<<type_crossover<< endl;

	for (int instance=0;instance<n_instances;instance++) {	
		cout <<"Instance: "<< instance << endl;
		srand(instance);												// random seed  (for instance)	
		randnseed = -n_run-1;											// seed for gaussian random number generator. OBS.: must be (initially) negative	
		// Creating the NK Model
		PHI_Mk = new list<int>[N_Mk]; 
		Mk_model(N_Mk, model_Mk, PHI_Mk);								// modify matrix M (when necessary) and create PHI_Mk

		// Creating the Mk Model	
		Mk *Mk_instance = new Mk(N_Mk, N_Mk, PHI_Mk, M, func_num, OShift);				// from class Mk (cMk.h)
		//Mk_instance->print();											// print graph and costs of the elements for the NK Landscape
		for (int i_run=0;i_run<n_runs_per_instance;i_run++){
			cout <<"Run:"<< n_run <<", random seed: " << n_instances+n_run << endl;
			srand( n_instances+n_run );									// random seed  (for run)
			de(Mk_instance, n_run);					    				// run DE
			n_run++;
		}

		delete [] PHI_Mk;
		delete Mk_instance;
	}
	
	file_output(N_Mk,model_Mk,total_runs);					// save data

	// Desallocation of vectors and matrices	
	desaloc_matrixd (File_recComp,total_runs);
	desaloc_matrixd (File_best_ind,total_runs);
	delete [] time_run;
	delete [] file_gen;
	delete [] file_best_fitness;
	delete [] file_best_fitness_gen;
	delete [] file_sucRate;
	delete [] file_impRate;	
	delete [] vsort_aux;
	// Desallocation of vectors and matrices: CEC test suit
	delete [] M;
	delete [] OShift;
	delete [] x_bound;
	delete [] y_cec2017;
	delete [] z_cec2017;

		
	return 0;
}
