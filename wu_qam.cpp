//wu's bpsk ber


#define PI 3.1415926535897932
#define ROOT_2 1.414213562

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime> 
#include <math.h>
#include <string>
using namespace std;
//-----------------------------------------------------------
double q_func(double x);
double erfc(double x);
//----------------------------------------------------------
int randomgen();  //random    +1¡¬s and -1¡¬s
double gaussgen(double mean,double stddev);  //random    Gaussian(mean, std)

//----------------------------------------------------------
void QAMBer_theory(int big_m, double beg_ebno, double end_ebno, int num_pts);
//----------------------------------------------------------





int main(){
	srand(time(NULL));   //srand
	QAMBer_theory(16,0,15,10000);   //qamBer_theory
	
	
	
	
	
//----------------start------------------------------------	
	fstream output_QAMBer_exp;
	output_QAMBer_exp.open("QAMBer_exp.txt",ios::out);
	
	
	
	
	
	
	double beg_ebno=0,
		end_ebno=15,
		ebno_db=0,
		pts=60,
		errnum=0,
		d=1,
		Eav=10,
		M=16,
		
		
		receive_r0=0,
		receive_r1=0,
		
		random_r0=0,
		random_r1=0,
		
	//	detect=0,
		simBer=0,
		sigma=0,
		N=10000000,
		
		decis_r0=0,
		decis_r1=0;
		
	double delta_ebno = (end_ebno - beg_ebno)/(double)(pts-1);
	double ebno_numeric=0;
	
	
	
	
	for(int i=0; i<pts; i++){
		ebno_db = beg_ebno + i * delta_ebno;
    	ebno_numeric = pow(10.0, ebno_db/10.0);
    	sigma=sqrt(Eav)/sqrt(8*ebno_numeric);
    	
    	for(int j=0;j<N;j++){
		
			random_r0=randomgen();  //   random    +3¡¬s , +1¡¬s , -1¡¬s , -3¡¬s
			random_r1=randomgen();
			
 			receive_r0=random_r0+gaussgen(0,sigma);   //add gauss
 			receive_r1=random_r1+gaussgen(0,sigma); 
 			
 			
 			if(receive_r0>2){
 				decis_r0=3;
			 }
			 else if(receive_r0>0){
			 	decis_r0=1;
			 }
 			 else if(receive_r0>-2){
			 	decis_r0=-1;
			 }
 			 else{
 			 	decis_r0=-3;
			}
 			
 			if(receive_r1>2){
 				decis_r1=3;
			 }
			 else if(receive_r1>0){
			 	decis_r1=1;
			 }
 			 else if(receive_r1>-2){
			 	decis_r1=-1;
			 }
 			 else{
 			 	decis_r1=-3;
			}
 			
 			
			
		//	if((decis_r0-receive_r0)==4||(decis_r0-receive_r0)==-4)  {errnum+=2;}
		    if((decis_r0!=random_r0)||(decis_r1!=random_r1))  {errnum++;}
			  
//			cout << random_r0 << ", " <<  receive_r0  <<", " <<  decis_r0  <<endl ;
		//	if((decis_r1-receive_r1)==4||(decis_r1-receive_r1)==-4)  {errnum+=2;}
		//    if()  {errnum++;}
			
		}
	//	cout<<gaussgen(0,sigma)<<endl;
		simBer=errnum/N;
		simBer=simBer/4;
		output_QAMBer_exp << ebno_db << ", " << simBer  << endl;
	errnum=0;
	}
	
	output_QAMBer_exp.close();
	return 0;
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







int randomgen(){
	
	int randbit=0;
	randbit=2*(rand()%4)-3;
	return randbit;
}

 
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






void QAMBer_theory(int big_m, double beg_ebno, double end_ebno, int num_pts)
{
  double ebno_numeric, ebno_db, bit_err, symb_err;
  double m_factor ;
  fstream output;
  output.open("QAMBer_theory.txt",ios::out);
  
  double delta_ebno = (end_ebno - beg_ebno)/(double)(num_pts-1);
  
  m_factor = double(log2(big_m));
//sin_arg = PI/double(big_m);
  
  
  for(int n=0; n<num_pts; n++)
    {
    ebno_db = beg_ebno + n * delta_ebno;
    ebno_numeric = pow(10.0, ebno_db/10.0);
    symb_err = 3 * q_func(sqrt(0.2 * m_factor * ebno_numeric) );
    bit_err = symb_err/m_factor;
    output << ebno_db << ", " << bit_err << ", " << symb_err << endl;
    
    }
    //output << "x2=BpskBer_exp(:,1);y2=BpskBer_exp(:,2);semilogy(x2,y2,'k');hold on" << endl;
    output.close();

}
