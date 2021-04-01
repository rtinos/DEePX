/******************************************************************************\
*							Subfunctions for 	CEC17 Test Function Suite      *
* Adapted from:      														   *
*  CEC17 Test Function Suite for Single Objective Optimization				   *
\******************************************************************************/


// Linear transformation of the variables
double sr_sfunc(double *x, int nx, double *Os, double *Mr, double sh_rate, int l, list<int> *PHI_Mk) {
	int j;
	double x_new=0.0, x_phi_l;
	
	for(list<int>::iterator ii = PHI_Mk[l].begin(); ii != PHI_Mk[l].end(); ii++){
		j=*ii;
		x_phi_l=x[j]-Os[j]; 				// shift
		x_phi_l*=sh_rate;					// shrink to the original search range
		x_new+=Mr[l*nx+j]*x_phi_l;			// linear tranformation
	}
	
	return x_new;
}


// Bent Cigar
double bent_cigar_sfunc(double *x, int nx, double *Os, double *Mr, int l, list<int> *PHI_Mk) {
    double f, x_new;
        
	x_new=sr_sfunc (x, nx, Os, Mr, 1.0, l, PHI_Mk); 	// transform x to x_new

	if (l==0)
		f = x_new*x_new;
	else
		f = pow(10.0,6.0)*x_new*x_new;
	
	return f;

}


// Rosenbrock's (observation: only shift  is applied here)
double rosenbrock_sfunc(double *x, int nx, double *Os, int l) 
{
	double tmp1, tmp2, f=0.0, sh_rate, x1, x2;
	
	if (l==nx-1)
		return f;
	
	sh_rate=2.048/100.0;
	
	x1=x[l];
	x1-=Os[l]; 				// shift
	x1*=sh_rate;				// shrink to the original search range
	x1+=1.0;					// shift to orgin
	x2=x[l+1];
	x2-=Os[l+1]; 				// shift
	x2*=sh_rate;				// shrink to the original search range
	x2+=1.0;					// shift to orgin

	tmp1=x1*x1-x2;
	tmp2=x1-1.0;
	f=100.0*tmp1*tmp1+tmp2*tmp2;

	return f;
}


// Rastrigin’s 
double rastrigin_sfunc(double *x, int nx, double *Os, double *Mr, int l, list<int> *PHI_Mk) {
    double f, x_new;
        
	x_new=sr_sfunc (x, nx, Os, Mr, 5.12/100.0, l, PHI_Mk); 	// transform x to x_new

	f = ( x_new*x_new - 10.0*cos(2.0*PI*x_new) + 10.0);
	
	return f;

}


// Schwefel’s 
double schwefel_sfunc(double *x, int nx, double *Os, double *Mr, int l, list<int> *PHI_Mk) {
    double f=0.0, x_new, tmp;
        
	x_new=sr_sfunc (x, nx, Os, Mr, 1000.0/100.0, l, PHI_Mk); 	// transform x to x_new

	x_new += 4.209687462275036e+002;
	if (x_new>500){
		f-=(500.0-fmod(x_new,500))*sin(pow(500.0-fmod(x_new,500),0.5));
		tmp=(x_new-500.0)/100;
		f+= tmp*tmp/nx;
	}
	else if (x_new<-500){
		f-=(-500.0+fmod(fabs(x_new),500))*sin(pow(500.0-fmod(fabs(x_new),500),0.5));
		tmp=(x_new+500.0)/100;
		f+= tmp*tmp/nx;
	}
	else
		f-=x_new*sin(pow(fabs(x_new),0.5));
	
	f+=4.189828872724338e+002;
	
	return f;

}


// compute evaluation for subfunctions
double cec17_test_sfunc( double *x, double *M, double *OShift, int nx, int f_num, int l, list<int> *PHI_Mk){
	double f;


	switch(func_num){
		case 1:	
			f=bent_cigar_sfunc(x,nx,OShift,M,l,PHI_Mk);
			break;
		case 4:	
			f=rosenbrock_sfunc(x,nx,OShift,l);
			break;
		case 5:
			f=rastrigin_sfunc(x,nx,OShift,M,l,PHI_Mk);
			break;
		case 10:	
			f=schwefel_sfunc(x,nx,OShift,M,l,PHI_Mk);
			break;
		default:
			cout<<"Incorrect argument!"<<endl;
			cout<<"Choices for CEC17 Test Function number: 1-Bent Cigar;";
			cout<<" 4-Rosenbrock's Function; ";
			cout<<" 5-Rastrigin's Function; ";
		    cout<<" 10-Schwefel's Function; ";	
			cout<<endl;
			exit(1);
			break;
	}
		
	return (-f);

}


// compute evaluation for subfunctions
double cec17_test_bias( int f_num ){
	double bias;

	switch(func_num){
		case 1:	
			bias=100.0;
			break;
		case 4:	
			bias=400.0;
			break;
		case 5:
			bias=500.0;
			break;
		case 10:	
			bias=1000.0;
			break;
		default:
			cout<<"Incorrect argument!"<<endl;
			cout<<"Choices for CEC17 Test Function number: 1-Bent Cigar;";
			cout<<" 4-Rosenbrock’s Function; ";
			cout<<" 5-Rastrigin’s Function; ";
		    cout<<" 10-Schwefel’s Function; ";	
			cout<<endl;
			exit(1);
			break;
	}		
	return bias;
	
}


// compute fitness (considering minimization)
double cec17_compFitness( double *x, double *M, double *OShift, int nx, int f_num, list<int> *PHI_Mk){
	double f=0.0;

	for (int l=0;l<nx;l++)
		f+=cec17_test_sfunc(x,M,OShift,nx,f_num,l,PHI_Mk);
	f-=cec17_test_bias(f_num);
	
	return f;

}

