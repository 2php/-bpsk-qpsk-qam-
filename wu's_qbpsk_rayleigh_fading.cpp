//wu's dbpsk rayleigh fading  
#define PI 3.1415926535897932
#define ROOT_2 1.414213562

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime> 
#include <cmath>
#include <string>
#include <complex>
#include <vector>
using namespace std;

double erfc(double x){
	double w, y, z;
	z = fabs(x);
	w = 1.0 / (1.0 + 0.5*z);
	y = w*exp(-z*z-1.26551223+w*(1.00002368+w*(0.37409196+
		w*(0.09678418+w*(-0.18628806+w*(0.27886807+w*(-1.13520398+
        w*(1.48851587+w*(-0.82215223+w*0.17087277)))))))));
    if(x<0.0) y = 2.0 - y;
    return(y);
}	

double q_func(double x){
	double value;
	if(x==0){
		value=0.5;
	}
	else{
		value=0.5*erfc(x/ROOT_2); 
	}
	return value;
}


//----------------------------------------------------------
int bpsk_randomgen(){
	
	int randbit=0;
		
	randbit = 2*(rand()%2)-1;
	return randbit;
} //random    +1¡¬s and -1¡¬s
	


 
double gaussgen(double mean, double stddev)
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
 //random    Gaussian(mean, std)


complex<double> h(double wm, double t ){

	double N0 = 8.0;
	double N=4*N0+2;
	double n=0.0;
	complex<double> h , hh ;
	
	hh.real()=0.0;
	hh.imag()=0.0;
	for(n=0;n<N0;n++){
		hh.real() += cos(wm*t*cos((2*PI*n/N)+(PI/2*N)));
		hh.imag() += sin(wm*t*sin((2*PI*n/N)+(PI/2*N)));
	}
	h.real()=(1/N0) * hh.real();
	h.imag()=(1/N0) * hh.imag();
							 
	return h;
	
}






int main(){
	
	fstream output_dqpsk_theory,output_dqpsk_exp;
	output_dqpsk_theory.open("output_dqpsk_theory.txt",ios::out);
	output_dqpsk_exp.open("output_dqpsk_exp.txt",ios::out);
	
	srand(time(NULL));
	double Ts=0.00003;
	double T=10000;
	
	unsigned long int  bpsk_length=25000;
	int bpsk_I[bpsk_length];
	int dpsk_I[bpsk_length+1];
	complex <double> Y1_I[bpsk_length+1] ;
	double  R1_I[bpsk_length];
	int detect_I[bpsk_length];
	
	int bpsk_Q[bpsk_length];
	int dpsk_Q[bpsk_length+1];
	complex <double> Y1_Q[bpsk_length+1] ;
	double  R1_Q[bpsk_length];
	int detect_Q[bpsk_length];
	
	
	double SNR_db = 0;
	double err=0;
	double sum_err=0;
	double SNR_exp;
	double sigma;
	
	
	for(SNR_db=0;SNR_db<20;SNR_db+=0.1){
		
		SNR_exp = pow(10,(SNR_db/10));
		sigma=1/sqrt(8*SNR_exp);
		
		
			
		
		for(int i=0;i<bpsk_length;i++){
			bpsk_I[i]= bpsk_randomgen(); 
			bpsk_Q[i]= bpsk_randomgen();
		}
		dpsk_I[0]=1;
		dpsk_Q[0]=1;
		
		for(int i=0;i<bpsk_length;i++){
			
	    	dpsk_I[i+1] = dpsk_I[i] * bpsk_I[i];
	    	dpsk_Q[i+1] = dpsk_Q[i] * bpsk_Q[i];
	    	
	    	Y1_I[i+1].real() = h(2*PI*80 , i*Ts).real()*dpsk_I[i+1]  + gaussgen (0,sigma )/1.414;
	    	Y1_I[i+1].imag() = h(2*PI*80 , i*Ts).imag()*dpsk_I[i+1]  + gaussgen (0,sigma )/1.414;
	    	
	    	Y1_Q[i+1].real() = h(2*PI*80 , i*Ts).real()*dpsk_Q[i+1] +  gaussgen (0,sigma )/1.414;
	    	Y1_Q[i+1].imag() = h(2*PI*80 , i*Ts).imag()*dpsk_Q[i+1]  + gaussgen (0,sigma )/1.414;


	    	
	 	}
	 	
	 	
	 	
	 	
	 	for(int i=0;i<bpsk_length;i++){
	 		
	 		
	 		
	    	R1_I[i] = Y1_I[i+1].real() * Y1_I[i].real() + Y1_I[i+1].imag() * Y1_I[i].imag();
	    	R1_Q[i] = Y1_Q[i+1].real() * Y1_Q[i].real() + Y1_Q[i+1].imag() * Y1_Q[i].imag();
	    	if(R1_I[i]>0) detect_I[i] = 1;
	    	else detect_I[i] = -1;
	    	
	    	if(R1_Q[i]>0) detect_Q[i] = 1;
	    	else detect_Q[i] = -1;
	//    	cout<<detect[i]<<endl;
	    	if((detect_I[i]!=bpsk_I[i])||(detect_Q[i]!=bpsk_Q[i]))  err++;
	    	
	    	
	    	
		 }
 		sum_err = (double)err/bpsk_length;
		sum_err/=2;
	 	output_dqpsk_exp << sum_err << ", " << SNR_db << endl;
	 	
	 	
 		err = 0;
 		sum_err=0;
 		
	}
 	








	
//-------------------------------------------
//-------------------------------------------
//-------------------------------------------
//-------------------------------------------
//-------------------------------------------


	
	double SNR=0.1 , snr123 = 0;
	double bessel_term = j0(2*PI*80*Ts);
	double oneover_SNR =1.0/SNR; 
 	for(int i=0;i<2000;i++){
	 	SNR += 0.01;
	 	snr123 = pow(10,(SNR/10));
	 	oneover_SNR=1.0/snr123;
//	 	cout<< oneover_SNR<<endl;
		output_dqpsk_theory << 0.5-0.5*bessel_term/sqrt(2+oneover_SNR+oneover_SNR-oneover_SNR*oneover_SNR*0.5-bessel_term)<< ", " << SNR <<endl;
	}

	output_dqpsk_exp.close();
	output_dqpsk_theory.close();
	return 0;  
}    
