/******************************************************************************\
*				  				 Files Manipulation							 *
\******************************************************************************/

#include "defs.h"
#include<cstring>
#include<fstream>


/******************************************************************************\
* 					Save data : end of the simulation						   *
\******************************************************************************/
void file_output(int N, int model_Mk, int total_runs)
{
	char *name_p;
	char name[CHAR_LEN];
	FILE *Bestfit_file, *Bestind_file, *Time_file, *Gen_file, *sucRate_file, *impRate_file;

    name_p = name;
		
  	// Best fitness in each generation for run 0
  	if (save_datagen_flag==1){
  		FILE *Bfg_file;
		sprintf(name,"bfg_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
		if ((Bfg_file = fopen(name_p,"w"))==NULL) {
			puts("The file bfg to be saved cannot be open \n");
			exit(1);
		}
		for (int i=0;i<max_gen;i++) {
			fprintf(Bfg_file,"%.5f ",file_best_fitness_gen[i]);
		}
		fclose(Bfg_file);
	}

    // Best fitness 
	sprintf(name,"bfi_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
	if ((Bestfit_file = fopen(name_p,"w"))==NULL) {
		puts("The file bfi to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Bestfit_file,"%.10f ",file_best_fitness[i]);
	}
	fclose(Bestfit_file);
		
	 // Best individuals
	sprintf(name,"bind_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
	if ((Bestind_file = fopen(name_p,"w"))==NULL) {
		puts("The file bind to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		for (int gene=0;gene<chrom_size;gene++)
			fprintf(Bestind_file,"%.5f ",File_best_ind[i][gene]);
		fprintf(Bestind_file,"\n");
	}
	fclose(Bestind_file);
	
  	// Time for each run
	sprintf(name,"time_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
	if ((Time_file = fopen(name_p,"w"))==NULL) {
		puts("The file time to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Time_file,"%.2f ",time_run[i]);
	}
	fclose(Time_file);
	
	// Number of generations for each run
	sprintf(name,"gen_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
	if ((Gen_file = fopen(name_p,"w"))==NULL) {
		puts("The file gen to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++) {
		fprintf(Gen_file,"%d ",file_gen[i]);
	}
	fclose(Gen_file);		
		
	// Crossover Statistics: successful recombination rate 
	sprintf(name,"sucRate_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
	if ((sucRate_file = fopen(name_p,"w"))==NULL) {
		puts("The file sucRate to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++)		
		fprintf(sucRate_file,"%.8f ",file_sucRate[i]);		
	fclose(sucRate_file);
	
	// Crossover Statistics: improvement rate 
	sprintf(name,"impRate_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
	if ((impRate_file = fopen(name_p,"w"))==NULL) {
		puts("The file impRate to be saved cannot be open \n");
		exit(1);
	}
	for (int i=0;i<total_runs;i++)		
		fprintf(impRate_file,"%.8f ",file_impRate[i]);		
	fclose(impRate_file);
	
	if (type_crossover==3 || type_crossover==4 || type_crossover==5){
			FILE *recComp_file;
			sprintf(name,"recComp_f%d_N%d_m%d_c%d.dat",func_num,N,model_Mk,type_crossover);
			if ((recComp_file = fopen(name_p,"w"))==NULL) {
				puts("The file recComp to be saved cannot be open \n");
				exit(1);
			}
			for (int i=0;i<total_runs;i++)	
				fprintf(recComp_file,"%.3f ",File_recComp[i][0]);
			fprintf(recComp_file,"\n");
			for (int i=0;i<total_runs;i++)	
				fprintf(recComp_file,"%.0f ",File_recComp[i][1]);				
			fclose(recComp_file);		
	}
	
}
