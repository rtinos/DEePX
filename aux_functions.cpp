/******************************************************************************\
*								Diverse Auxiliary Functions					   *
\******************************************************************************/
#include "defs.h"
#include <cstdlib>
#include <cmath>
 
/************************************************************************************************\
*	Routines book: Numerical Recipes in C														*
\************************************************************************************************/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/**********************************************************************************************\
*	 Minimal ”random num er generator of Park and Miller with Bays-Durham shu .e and added      *
*  safeguards.Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint   *
*  values).Call with idum a negative integer to initialize;thereafter,do not alter idum between *
*  successive deviates in a sequence. RNMX should approximate the largest floating value that is *
*  less than 1. 																										  *
\********************************************************************************************/
double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy) { 													//Initialize.
		if (-(*idum) < 1) *idum=1; 												//Be sure to prevent idum =0
   	else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { 													//Load the shu .e ta le (after 8 warm-ups).
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
   		if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
		}
	k=(*idum)/IQ; 																	//Start here when not initializing.
	*idum=IA*(*idum-k*IQ)-IR*k; 												// Compute idum=(IA*idum) % IM without over-
																					// .ows y Schrage ’s method.
   if (*idum < 0) *idum += IM;
   j=iy/NDIV; 							//Will e in the range 0..NTAB-1
   iy=iv[j]; 							//Output previously stored value and re .ll the
   iv[j] = *idum; 					//shu .e ta le.
   if ((temp=AM*iy) > RNMX) return RNMX; //Because users don ’t expect endpoint values.
   else return temp;
}

/*******************************************************************************\
*	 Returns a normally distributed deviate with zero mean and unit variance,    *
*	using ran1(idum) as the source of uniform deviates. (Numerical recipes in C) *
\******************************************************************************/
double gasdev(long int *idum)
{
	double ran1(long int *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

   if (*idum < 0) iset=0; 												//Reinitialize.
	if (iset == 0) { 														//We don’t have an extra deviate handy, so
		do {
			v1=2.0*ran1(idum)-1.0; 												//pick two uniform numbers in the square ex-
			v2=2.0*ran1(idum)-1.0;												//tending from -1 to +1 in each direction,
 			rsq=v1*v1+v2*v2;														// see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); 								// and if they are not, try again.
		fac=sqrt(-2.0*log(rsq)/rsq);
		// Now make the Box-Muller transformation to get two normal deviates. Return one and
		//save the other for next time.
		gset=v1*fac;
		iset=1; 																//Set flag.
		return v2*fac;
	} else { 																//We have an extra deviate handy,
   	iset=0; 																//so unset the flag,
		return gset; 														  //and return it.
	}
}


/******************************************************************************\
*								 Dynamic Allocation: Matrix of Integers		   *
\******************************************************************************/
int **aloc_matrixi(int lines , int collums)
{
	int i, **Matrix;
	
	Matrix = new int*[lines];
	for (i=0;i<lines;i++) {
		Matrix[i] = new int[collums];
	}
	if (!Matrix) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}

	return Matrix;
}


/******************************************************************************\
*								 Dynamic Allocation: Matrix of Doubles		   *
\******************************************************************************/
double **aloc_matrixd(int lines , int collums)
{
	double **Matrix;
	
	Matrix = new double*[lines];
	for (int i=0;i<lines;i++) {
		Matrix[i] = new double[collums];
	}
	if (!Matrix) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}

	return Matrix;
}


/******************************************************************************\
*								Dynamic Allocation: Vector of Integers		   *
\******************************************************************************/
int *aloc_vectori(int lines)
{
	int *vector;

	vector = new int[lines];
	if (!vector) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}
	return vector;
}


/******************************************************************************\
*								Dynamic Allocation: Vector of Doubles		   *
\******************************************************************************/
double *aloc_vectord(int lines)
{
	double *vector;

	vector = new double[lines];
	if (!vector) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}
	return vector;
}


/******************************************************************************\
*								Dynamic Allocation: Vector of individuals						 *
\******************************************************************************/
individual *aloc_vectorind(int lines)
{
	individual *vector;

	vector = new individual[lines];
	if (!vector) {
		cout<<"Allocation Error!"<<endl;
		exit(1);
	}
	return vector;
}


/******************************************************************************\
*								 Dynamic Desallocation: Matrix of Integers	   *
\******************************************************************************/
void desaloc_matrixi(int **Matrix , int lines)
{
	int i;

	for(i=0;i<lines;i++) {
		delete [] Matrix[i];
	}
	delete [] Matrix;

}


/******************************************************************************\
*								 Dynamic Desallocation: Matrix of Doubles	   *
\******************************************************************************/
void desaloc_matrixd(double **Matrix , int lines)
{

	for(int i=0;i<lines;i++) {
		delete [] Matrix[i];
	}
	delete [] Matrix;

}


/******************************************************************************\
*						 Random Integer between L_range and H_range			   *
\******************************************************************************/
int random_int(int L_range, int H_range)
{
	return(  (int) ( (rand()/(RAND_MAX+1.0))*(H_range-L_range+1)+L_range ) );  	// random integer beteween [L_range and H_range]
}


/******************************************************************************\
*								 Random double in [0.0,1.0]			 		   *
\******************************************************************************/
double random_dou(void)
{
	return(  rand() / double(RAND_MAX) );  //  random double in [0.0, 1.0]:
}



/******************************************************************************\
*		 Multiplication of matrices M=AB : double						 	   *
\******************************************************************************/
void multMatrix(double **M, double **A, int l_A, int c_A, double **B, int l_B, int c_B)
{
	int i ,j, k;

	if (c_A!=l_B) {
		cout<<"Error - multiplication: size of the matrices!"<<endl;
		exit (1);
	}

	for (i=0;i<l_A;i++) {
		for (j=0;j<c_B;j++) {
			M[i][j] = 0.0;
			for (k=0;k<c_A;k++) {
					M[i][j] = M[i][j]+ A[i][k]*B[k][j];
			}
		}
	}
}


/******************************************************************************\
*		Transpose of matrix M : double								 		   *
\******************************************************************************/
void transpose(double **Mt, double **M , int l , int c)
{
	int i , j;

	for (i=0;i<l;i++) {
		for (j=0;j<c;j++) {
			Mt[j][i]=M[i][j];
		}
	}
}


/******************************************************************************\
* check if doubles are equal												   *
\******************************************************************************/
int isEqual(double x, double y)
{
 
	if ( fabs(x-y)>EPS1 )
		return (0);
	else
		return (1);
	
}


/******************************************************************************\
*		 print matrix ( double)						 	   					   *
\******************************************************************************/
void printMatrix(double **A, int l, int c){
	
	cout<<"Matrix:"<<endl;
	cout<<fixed;
	for (int i=0;i<l;i++) {
		for (int j=0;j<l;j++) 
			cout<<A[i][j]<<"\t";
		cout<<endl;
	}
	
}


/******************************************************************************\
*		 rotation matrix ( rotation - adjacent with sobreposition)			   *
\******************************************************************************/
void rotMatrix1(double **A, int l)
{
	int i, j;
	double rho=0.5, theta, cos_theta, sin_theta, **R1, **R2;
	
	if ((l%2)!=0) {
		cout<<"Error - rotMatrix (l must be even )!"<<endl;
		exit (1);
	}
	R1=aloc_matrixd(l,l);
	R2=aloc_matrixd(l,l);
	
	for (i=0;i<l;i++) {
		R1[i][i]=1.0;
		R2[i][i]=1.0;
		for (j=i+1;j<l;j++) {
			R1[i][j]=0.0;
			R1[j][i]=0.0;
			R2[i][j]=0.0;
			R2[j][i]=0.0;		
		}		
	}
		
	// rotation matrix R1
	for (i=0;i<l-1;i=i+2) {
		j=i+1;
		theta = random_dou()*rho*PI;
		cos_theta = cos(theta);
		sin_theta = sin(theta);
		R1[i][i] = cos_theta;
		R1[i][j] = -sin_theta;
		R1[j][i] = sin_theta;
		R1[j][j] = cos_theta;
	}	
			
	// rotation matrix R2
	for (i=1;i<l-2;i=i+2) {
		j=i+1;
		theta = random_dou()*rho*PI;
		cos_theta = cos(theta);
		sin_theta = sin(theta);
		R2[i][i] = cos_theta;
		R2[i][j] = -sin_theta;
		R2[j][i] = sin_theta;
		R2[j][j] = cos_theta;
	}
	theta = random_dou()*rho*PI;
	cos_theta = cos(theta);
	sin_theta = sin(theta);
	R2[l-1][l-1]=cos_theta;
	R2[l-1][0]=-sin_theta;
	R2[0][l-1]=sin_theta;
	R2[0][0]=cos_theta;		
	
	// A=R1*R2
	multMatrix(A, R1, l, l, R2, l, l);
	// printMatrix(A, l, l);
	
	// Checking if it is orthogonal and epistasis degree
	/*	double **Ort, **At;
		int flag_ort=1, flag_k=1, k_aux;
		Ort=aloc_matrixd(l,l);
		At=aloc_matrixd(l,l);
		transpose(At, A , l, l);
		multMatrix(Ort, A, l, l, At, l, l);	
		for (i=0;i<l;i++){
			k_aux=0;
			for (j=0;j<l;j++){
				if (isEqual(A[i][j],0.0)==0)
					k_aux++;
				if ( i!=j){
				 	if (isEqual(Ort[i][j],0.0)==0)
						flag_ort=0;
				}
			    else{
					if (isEqual(Ort[i][j],1.0)==0)
			    		flag_ort=0;
			    }
			}
			if (k_aux!=4)
				flag_k=0;		
		}
		desaloc_matrixd(At,l);
		desaloc_matrixd(Ort,l);
		if (flag_k==0 || flag_ort==0){
			cout<<"Error - rotMatrix (not orthogonal or k!=4)!"<<endl;
			exit (1);
		}
	*/
				
	desaloc_matrixd(R1,l);
	desaloc_matrixd(R2,l);
}


/******************************************************************************\
* Random Permutation of a vector of integers (only size first elements of inp) *
\******************************************************************************/
void rand_perm_size(int *inp, int *out, int size_inp, int size_out)
{
	int j, aux, *auxv;
	
	auxv=aloc_vectori(size_inp);
	
	for(int i=0;i<size_inp;i++)
		auxv[i]=inp[i];
	
	for(int i=0;i<size_out;i++) {
		j= random_int(i,size_inp-1);  		
		aux=auxv[i];
		auxv[i]=auxv[j];
		auxv[j]=aux;
		out[i]=auxv[i];
	}
	
	delete [] auxv;
}

/******************************************************************************\
* check if integer vectors are equal 											*
\******************************************************************************/
int isVectorEqual(int *v1, int *v2, int size_v)
{
	int i=0;

	while (i<size_v && v1[i]==v2[i])
		i++;
	
	if (i<size_v)
		return (0);
	else
		return (1);
	
}


