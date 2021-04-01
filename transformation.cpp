#include "defs.h"

/*******************************************************************************\
*	 		Exponential Crossover (similar to 2-point Crossover	in GAs) 	    *
*	ref: Zaharie, D. (2009). Influence of crossover on the behavior of 	 	    *
*    differential evolution algorithms. Applied soft computing, 9(3), 1126-1138.*	 						 *
\*******************************************************************************/
void Point2X(allele *parent1, allele *parent2, allele *offspring){
		int p1, p2, aux;
		double aux_d;

		// defining the crossover points
		p1 =random_int(0,chrom_size-1);		// point 1
		p2 =random_int(0,chrom_size-1);		// point 2
		while (p1 == p2)
			p2 =random_int(0,chrom_size-1);	// point 2
		if (p1>p2) {
			aux=p1;
			p1=p2;
			p2=aux;
		}											 
		// generating the offspring
		aux_d=random_dou();
		if (aux_d<0.5){
			for (int gene=0;gene<chrom_size;gene++) {
				if (gene<p1 || gene>=p2) 
					offspring[gene] = parent1[gene];	
				else
					offspring[gene] = parent2[gene];
			}
		}
		else{
			for (int gene=0;gene<chrom_size;gene++) {
				if (gene<p1 || gene>=p2) 
					offspring[gene] = parent2[gene];	
				else
					offspring[gene] = parent1[gene];
			}
		}
		
}


/*******************************************************************************\
* 			Binomial Crossover (similar to Uniform Crossover in GAs)	   	    *
*	ref: Zaharie, D. (2009). Influence of crossover on the behavior of 	 	    *
*    differential evolution algorithms. Applied soft computing, 9(3), 1126-1138.*	 						 *
\*******************************************************************************/
void UX(allele *parent1, allele *parent2, allele *offspring  ){
		int jrand;
		double aux, CR_DE=0.9;
		
		jrand=random_int(0,chrom_size-1);
		for (int gene=0;gene<chrom_size;gene++){
			aux=random_dou();							
			if (aux<CR_DE || gene==jrand)
				offspring[gene] = parent1[gene];			
			else
				offspring[gene] = parent2[gene];
		}	
				
}




