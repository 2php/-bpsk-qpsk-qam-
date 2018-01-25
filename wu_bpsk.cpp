//wu's bpsk ber


#define PI 3.1415926535897932
#define ROOT_2 1.414213562

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime> 
#include <cmath>
#include <string>
#include <iomanip> 

using namespace std;
//-----------------------------------------------------------
long double q_func(long double x);
long double erfc(long double x);
//----------------------------------------------------------
int randomgen();  //random    +1¡¬s and -1¡¬s
//double rand_normal(double mean,double std);  //random    Gaussian(mean, std)
double rand_normal(double mean, double std);
//----------------------------------------------------------
void BpskBer_theory( double beg_ebno,  double end_ebno, int num_pts);
//----------------------------------------------------------





int main(){
	srand(time(NULL));   //srand
	BpskBer_theory(0,15,10000);   //BpskBer_theory
	
	
	
	
	
//----------------start------------------------------------	
	fstream output_BpskBer_exp;
	output_BpskBer_exp.open("BpskBer_exp.txt",ios::out);
	
	long double beg_ebno=0,
		end_ebno=15,
		ebno_db=0,
		pts=60, 
		errnum=0,
		receive=0,
		
		simBer=0,
		sigma=0,
		N=10000000;
	long double delta_ebno = (end_ebno - beg_ebno)/(double)(pts-1);
	long double ebno_numeric=0;
	int random=0,iphat = 0;
	int Eb = 1;
	
	
	
	for(int i=0; i<pts; i++){
		ebno_db = beg_ebno + i * delta_ebno;
    	ebno_numeric = pow(10.0, ebno_db/10.0);
    	sigma=Eb/sqrt(2*ebno_numeric);
    	
    	for(int j=0;j<N;j++){
		
			random=randomgen();  //   random    +1¡¬s and -1¡¬s
 			receive=random + rand_normal(0,sigma);   //add gauss
 			if(receive>0) iphat = 1;
 			else iphat = -1;
			if(iphat ^ random) errnum++;      //detection and counting the number of errors
		}
			
		simBer=errnum/N;
		cout <<  errnum<< endl;
		output_BpskBer_exp << ebno_db << ", " << simBer << endl;
	errnum=0;
	}

	output_BpskBer_exp.close();
	return 0;
}












long double q_func(long double x){
	long double value;
	if(x==0){
		value=0.5;
	}
	else{
		value=0.5*erfc(x/ROOT_2); 
	}
	return value;
}






long double erfc(long double x){
	long double w, y, z;
	z = fabs(x);
	w = 1.0 / (1.0 + 0.5*z);
	y = w*exp(-z*z-1.26551223+w*(1.00002368+w*(0.37409196+
		w*(0.09678418+w*(-0.18628806+w*(0.27886807+w*(-1.13520398+
        w*(1.48851587+w*(-0.82215223+w*0.17087277)))))))));
    if(x<0.0) y = 2.0 - y;
    return(y);
}	







int randomgen(){
	
	int randbit=0;
	randbit=2*(rand()%2)-1;
	return randbit;
}
 /*
 double rand_normal(double mean,double std){


	double  u;
	double v;
	double r;
	do{
		u = rand() / (double)RAND_MAX;
		v = rand() / (double)RAND_MAX;
		//r = u*u + v*v;
	}while(u  == 0.0|| u == 1);
	
    
    double x = sqrt(-2 * log(u)) * cos(2 * PI * v) * std + mean;
    
	return x;
}

 */
 
 
double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}






 
void BpskBer_theory(double beg_ebno, double end_ebno, int num_pts)
{
  double ebno_numeric, ebno_db, bit_err;
  fstream output;
  output.open("BpskBer_theory.txt",ios::out);
  
  double delta_ebno = (end_ebno - beg_ebno)/(double)(num_pts-1);
  for(int n=0; n<num_pts; n++)
    {
    ebno_db = beg_ebno + n * delta_ebno;
    ebno_numeric = pow(10.0, ebno_db/10.0);
    bit_err = q_func(sqrt(2.0*ebno_numeric));
    output << ebno_db << ", " << bit_err << endl;
    
    }
  //  output << "x2=BpskBer_exp(:,1);y2=BpskBer_exp(:,2);semilogy(x2,y2,'k');hold on" << endl;
    output.close();
 
}
