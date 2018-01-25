//wu's qpsk ber


#define PI 3.1415926535897932
#define ROOT_2 1.414213562

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime> 
#include <math.h>
#include <string>
#include <complex>
using namespace std;
//-----------------------------------------------------------
double q_func(double x);
double erfc(double x);
//----------------------------------------------------------
int randomgen();  //random    +1¡¬s and -1¡¬s
double gaussgen(double mean,double stddev);  //random    Gaussian(mean, std)

//----------------------------------------------------------
void QpskBer_theory(double beg_ebno, double end_ebno, int num_pts);
//----------------------------------------------------------





int main(){
	srand(time(NULL));   //srand
 	QpskBer_theory(0,15,10000);   //BpskBer_theory
	
	
	
	
	
//----------------start------------------------------------	
	fstream output_QpskBer_exp;
	output_QpskBer_exp.open("QpskBer_exp.txt",ios::out);
	
	double beg_ebno=0,
		end_ebno=15,
		ebno_db=0,
		pts=60,
		errnum=0,
		
		receive_r0=0,
		receive_r1=0,
		
		randominput=0,
		decis=0,
		
		Ber=0,
		sigma=0,
		N=10000000;
	double delta_ebno = (end_ebno - beg_ebno)/(double)(pts-1);
	double ebno_numeric=0;
	
	
	
	
	for(int i=0; i<pts; i++){
		ebno_db = beg_ebno + i * delta_ebno;
    	ebno_numeric = pow(10.0, ebno_db/10.0);
    	sigma=1/sqrt(4*ebno_numeric);
    	
    	for(int j=0;j<N;j++){
		
			randominput=randomgen();  //   random    (+1,0)-->0 , (-1,0)-->1 , (0,+1)-->2 , (0,-1)-->3
			
			// matched filter outputs
			if(randominput==0){
				receive_r0= 1+gaussgen(0,sigma);  
				receive_r1= gaussgen(0,sigma);  
			}
			
			if(randominput==1){
				receive_r0= -1+gaussgen(0,sigma);  
				receive_r1= gaussgen(0,sigma);  
			}
 			
 			if(randominput==2){
				receive_r0= gaussgen(0,sigma);  
				receive_r1= 1+gaussgen(0,sigma);  
			}
			
			if(randominput==3){
				receive_r0= gaussgen(0,sigma);  
				receive_r1= -1+gaussgen(0,sigma);  
			}
			
			
			//detector
			if(receive_r0>abs(receive_r1))         {decis=0;}
			else if(receive_r0<-abs(receive_r1))   {decis=1;}
			else if(receive_r1>abs(receive_r0))    {decis=2;}
			else                                   {decis=3;}
			
			
			if(randominput!=decis)                  errnum++;
		}
		
		Ber=errnum/N;
		Ber /= 2; 
		output_QpskBer_exp << ebno_db << ", " << Ber << endl;
	errnum=0;
	}
	
	output_QpskBer_exp.close();
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
	randbit=(rand()%4);
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








void QpskBer_theory(double beg_ebno, double end_ebno, int num_pts)
{
  double ebno_numeric, ebno_db, bit_err,symb_err;
  fstream output;
  output.open("QpskBer_theory.txt",ios::out);
  
  double delta_ebno = (end_ebno - beg_ebno)/(double)(num_pts-1);
  for(int n=0; n<num_pts; n++)
    {
    ebno_db = beg_ebno + n * delta_ebno;
    ebno_numeric = pow(10.0, ebno_db/10.0);
    bit_err = q_func(sqrt(2.0*ebno_numeric));
    symb_err = bit_err * (2.0 - bit_err);
    output << ebno_db << ", " << bit_err << ", " << symb_err << endl;
    
    }
    //output << "x2=BpskBer_exp(:,1);y2=BpskBer_exp(:,2);semilogy(x2,y2,'k');hold on" << endl;
    output.close();

}

