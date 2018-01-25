//wu's enhanced modle rayleigh fading  
//one fader
#define PI 3.1415926535897932


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime> 
#include <cmath>
#include <string>
#include <complex>
#include <vector>
using namespace std;

int alpha_randomgen(){
	
	int ra=0;
	ra=rand()%2;
	return ra;
}


complex<double> h(double wm, double t ){
	//wn = max doppler spread
	//Ts = smapling period 
	//N = sample generate
	// one channel path
	double N0 = 8.0;
	double N=4*N0+2;
	double n=0.0;
//	double wn = wm * cos(2*PI*n/N);
//	double beta_n = n*PI/N;
//	double alpha_nk =  (2*PI*n/N)+(PI/2*N);
//	double phi_nk = n*PI/N0;
//	double phi_nk2 = PI*(n/N0-0.5-(1/2*N0));
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
	
//	fstream output_rayleigh_exp;
	fstream output_rayleigh_exp_distribution, output_rayleigh_exp_distribution_theory;
//	output_rayleigh_exp.open("rayleigh_exp.txt",ios::out);
	output_rayleigh_exp_distribution.open("rayleigh_exp_distribution.txt",ios::out);
	output_rayleigh_exp_distribution_theory.open("output_rayleigh_exp_distribution_theory.txt",ios::out);
	srand(time(NULL));
	double Ts=0.00003;
	double T=10000;
 
	for(double n=0 ;n<T;n++) {
	complex<double> output_h = h(2*PI*200,n*Ts);
	double atan_h =  atan2(output_h.imag(), output_h.real());
//	output_rayleigh_exp << abs(output_h) << ", " << atan_h  << ", "  <<  output_h.real() << ", "  <<  output_h.imag()  <<  ", "  <<  n*Ts    <<endl;
}



//------------------------------------------------
	
//	double tt=0.01;
	double begin=-1;
	double end=1;
//	double center_tt = begin + tt;;
	int index=0;
	double N = 10000000;
	int length=2000;
	double delta = (end-begin)/length;
	struct pdf{
		double data1;
		int data2;
	};
	struct pdf v[length];
	for(int i = 0;i<length ; i++){
		v[i].data2 = 0;
		v[i].data1 = 0;
	}
	   
//	for(index=0;index<length;index++){
		for(int i=0;i<N;i++){
			double num = h(2*PI*200,i*Ts).imag();
			index = ((int)((num-begin)*1000) );
	
			v[index].data2++;
			v[index].data1 = num;
	//		cout << index << endl;
		}
		for(int i=0;i<length;i++){
		output_rayleigh_exp_distribution << v[i].data1 << ", "<< v[i].data2/(N*delta) <<endl;

		}

//--------------------------------------------------------------
	
	
//	output_rayleigh_exp.close();
	output_rayleigh_exp_distribution.close();
	output_rayleigh_exp_distribution_theory.close();
	return 0;  
}    
