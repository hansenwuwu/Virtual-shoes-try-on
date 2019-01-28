#include "distribution.h"
#include "Statistic.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#define ERR 0.000001
#define DELTA 0.001

using namespace distribution;

//
// inline Functions
//

// factorial functions
inline unsigned long double fact(int x){
	if(x==0) return 1;
	else{
		unsigned long double sum;
		for(sum=1.0f; x!=0; x--) sum*=x; 
		return sum;
	}
}
// gamma function
inline unsigned long double gamma(double r){
	if(r==(int)r) return fact((int)r-1);
	double z=r-(int)r+1.0f;
	double prev=0.0f;
	double x=DELTA;
	unsigned long double output=4.0f*DELTA*pow(x, z-1.0f)*exp(-1.0f*x)/3;

	for(x+=DELTA; output-prev>ERR || x<1.0f; x+=DELTA*2){
		prev=output;
		output+=2.0f*DELTA*pow(x, z-1.0f)*exp(-1.0f*x)/3;
		output+=4.0f*DELTA*pow(x+DELTA, z-1.0f)*exp(-1.0f*(x+DELTA))/3;
	}
	output+=DELTA*pow(x, z-1.0f)*exp(-1.0f*x)/3;
	if(r<1){
		while(z>r){
			z--;
			output/=z;
		}
	}
	else{
		while(z<r){
			output*=z;
			z++;
		}
	}
	return output;
}
// Beta function: gamma(a)gamma(b)/gamma(a+b)
inline unsigned long double Beta(double a, double b){
	unsigned long double output=gamma(a);
	output/=gamma(a+b);
	return output*gamma(b);
}

//
// distribution class
//

DISTRIBUTION::DISTRIBUTION(void){
	_a=1.0f; _b=1.0f;
}
DISTRIBUTION::DISTRIBUTION(double a){
	_a=a;
}
DISTRIBUTION::DISTRIBUTION(double a, double b){
	_a=a; _b=b;
}
void DISTRIBUTION::set_a(double a){
	_a=a;
}
void DISTRIBUTION::set_b(double b){
	_b=b;
}

//
// Discrete Random Variables
//

// Uniform

Uniform_Disc::Uniform_Disc(void): DISTRIBUTION(0, 1){
}
Uniform_Disc::Uniform_Disc(double a, double b): DISTRIBUTION((int)a, (int)b){
	if(_b<_a) swap(_a, _b);
}
double Uniform_Disc::mean(void){
	return (_b+_a)/2;
}
double Uniform_Disc::var(void){
	return (SQUARE(_b-_a+1.0f)-1.0f)/12.0f;
}
double Uniform_Disc::sd(void){
	return sqrt(var());
}
double Uniform_Disc::f(int x){
	if(x>=_a && x<=_b) return 1.0f/(1+_b-_a);
	else return 0.0f;
}
double Uniform_Disc::f(double x){
	return f((int)x);
}
double Uniform_Disc::F(int x){
	if(x<_a) return 0.0f;
	else if(x>_b) return 1.0f;
	else return (x-_a+1)/(_b-_a+1);
}
double Uniform_Disc::F(double x){
	return F((int)x);
}
double Uniform_Disc::x(void){
	return _a+rand()%((int)_b-(int)_a+1);
}

// Poisson

Poisson::Poisson(void): DISTRIBUTION(1){
}
Poisson::Poisson(double ave): DISTRIBUTION(ave){
	if(ave<0.0f) {
		cerr<< "warning: parameter of Poisson distribution <0\n";
		_a=-_a;
	}
}
double Poisson::mean(void){
	return _a;
}
double Poisson::var(void){
	return _a*_a;
}
double Poisson::sd(void){
	return _a;
}
double Poisson::f(int x){
	if(x<0) return 0;
	double output=exp(-1.0f*_a);
	for(int i=1; i<=x; i++) {
		output*=_a;
		output/=(double)i;
	}
	return output;
}
double Poisson::f(double x){
	return f((int)x);
}
double Poisson::F(int x){
	if(x<0) return 0;
	double e=exp(-1.0f*_a);
	double sum, temp;
	int i, j;
	for(i=0, sum=0.0f; i<=x; i++, sum+=temp){
		for(j=1, temp=1; j<=i; j++) {
			temp*=_a;
			temp/=(double)j;
		}
	}
	return e*sum;
}
double Poisson::F(double x){
	return F((int)x);
}
double Poisson::x(void){
	double cdf=(double)rand()/RAND_MAX;
	double Fx=0.0f;
	int x=-1;
	while(Fx<cdf){
		x++;
		Fx+=f(x);
	}
	return (double)x;
}


//
// Continuous Random Variables
//

// Uniform

Uniform_Cont::Uniform_Cont(void): DISTRIBUTION(0.0f, 1.0f){
}
Uniform_Cont::Uniform_Cont(double a, double b): DISTRIBUTION(a, b){
	if(_b<_a) swap(_a, _b);
}
double Uniform_Cont::mean(void){
	return (_b+_a)/2;
}
double Uniform_Cont::var(void){
	return (SQUARE(_b-_a))/12.0f;
}
double Uniform_Cont::sd(void){
	return sqrt(var());
}
double Uniform_Cont::f(int x){
	return f((double)x);
}
double Uniform_Cont::f(double x){
	if(x>=_a && x<=_b) return 1.0f/(_b-_a);
	else return 0.0f;
}
double Uniform_Cont::F(int x){
	return F((double)x);
}
double Uniform_Cont::F(double x){
	if(x<_a) return 0.0f;
	else if(x>_b) return 1.0f;
	else return (x-_a)/(_b-_a);
}
double Uniform_Cont::x(void){
	return _a+(_b-_a)*(double)rand()/RAND_MAX;
}

// Normal Distribution

Normal::Normal(void): DISTRIBUTION(0.0f, 1.0f){
}
Normal::Normal(double ave, double var): DISTRIBUTION(ave, sqrt(fabs(var))){
}
double Normal::mean(void){
	return _a;
}
double Normal::var(void){
	return _b*_b;
}
double Normal::sd(void){
	return _b;
}
double Normal::f(int x){
	return f((double)x);
}
double Normal::f(double x){
	return 1.0f/sqrt(2.0f*M_PI)/_b*exp(-0.5f*SQUARE(x-_a)/SQUARE(_b));
}
double Normal::F(int x){
	return F((double)x);
}
double Normal::F(double x){
	double z=(x-_a)/_b;
	double F=0.5, i;
	if(z<-5.0) return 0.0;
	if(z>5.0) return 1.0;
	if(z==0.0) return 0.5;
	else if(z>0.0){
		for(i=0.0; i<z; i+=0.01)
			F+=0.005/sqrt(2.0*M_PI)*(exp(-0.5*SQUARE(i))+exp(-0.5*SQUARE(i+0.01)));
	}
	else{
		for(i=0.0; i>z; i-=0.01)
			F-=0.005/sqrt(2.0*M_PI)*(exp(-0.5*SQUARE(i))+exp(-0.5*SQUARE(i-0.01)));
	}
	return F+(z-i)/sqrt(2.0*M_PI)*exp(-0.5*SQUARE(z));
}
double Normal::x(void){
	// Box Muller Method
	double U1=(double)rand()/RAND_MAX;
	double U2=(double)rand()/RAND_MAX;
	double z=sqrt(-2.0*log(U1))*cos(2.0*M_PI*U2);
	return z*_b+_a;
}

// Exponential

Exponential::Exponential(void): DISTRIBUTION(1.0f){
}
Exponential::Exponential(double ave): DISTRIBUTION(1.0f/ave){
	if(ave<0.0f) {
		cerr<< "warning: parameter of exponential distribution <0.\n";
		_a=-_a;
	}
}
double Exponential::mean(void){
	return 1.0f/_a;
}
double Exponential::var(void){
	return 1.0f/_a/_a;
}
double Exponential::sd(void){
	return 1.0f/_a;
}
double Exponential::f(int x){
	return f((double)x);
}
double Exponential::f(double x){
	if(x<=0) return 0.0f;
	return _a*exp(-1.0f*_a*x);
}
double Exponential::F(int x){
	return F((double)x);
}
double Exponential::F(double x){
	if(x<=0) return 0.0f;
	return 1.0f-exp(-1.0f*_a*x);
}
double Exponential::x(void){
	double F=(double)rand()/RAND_MAX;
	return -1.0f*log(1.0f-F)/_a;
}

// Gamma Distribution

Gamma::Gamma(void): DISTRIBUTION(1.0f, 1.0f){
	_gamma=gamma(_b);
}
Gamma::Gamma(double scale, double shape): DISTRIBUTION(scale, shape){
	if(scale<0) _a=-_a;
	if(shape<0) _b=-_b;
	if(shape<0 || scale<0) cerr<< "warning: parameter of Gamma distribution <0.";
	_gamma=gamma(_b);
}
void Gamma::set_b(double b){
	_b=b;
	if(_b<0) _b=-_b;
	_gamma=gamma(_b);
}
double Gamma::mean(void){
	return _b/_a;
}
double Gamma::var(void){
	return _b/_a/_a;
}
double Gamma::sd(void){
	return sqrt(var());
}
double Gamma::f(int x){
	return f((double)x);
}
double Gamma::f(double x){
	if(x<=0.0f) return 0.0f;
	if((int)_b!=_b) return pow(_a, _b)*pow(x, _b-1.0f)*exp(-1.0f*_a*x)/_gamma;
	else{
		int i;
		long double sum=1.0f;
		for(i=0; i<_b; i++) sum*=_a;
		for(i=0; i<_b-1; i++) sum*=x;
		return sum*exp(-1.0f*_a*x)/_gamma;
	}
}
double Gamma::F(int x){
	return F((double)x);
}
double Gamma::F(double x){
	if(x<=0.0f) return 0.0f;
	double delta=x/5000;
	double F=4.0f*delta*f(delta)/3, k;
	for(k=2.0f*delta; k+2.0f*DELTA<x; k+=delta*2.0f){
		F+=delta*2.0f*f(k)/3;
		F+=delta*4.0f*f(k+delta)/3;
	}
	return F+(x-k)*f(x)/3;
}
double Gamma::x(void){
	double F=(double)rand()/RAND_MAX;
	double Fx=0.0f, x;
	for(x=DELTA; Fx<F; x+=DELTA){
		Fx+=DELTA*f(x);
	}
	return x-(Fx-F)*f(F);
}

// Erlang Distribution

Erlang::Erlang(void): DISTRIBUTION(1.0f, 1.0f){
	_gamma=gamma(_b);
}
Erlang::Erlang(double scale, double shape): DISTRIBUTION(scale, (int)shape){
	if(scale<0) _a=-_a;
	if(shape<0) _b=-_b;
	if(shape<0 || scale<0) cerr<< "warning: parameter of Erlang distribution <0.";
	_gamma=gamma(_b);
}
void Erlang::set_b(double b){
	_b=b;
	if(_b<0) _b=-_b;
	_gamma=gamma(_b);
}
double Erlang::mean(void){
	return _b/_a;
}
double Erlang::var(void){
	return _b/_a/_a;
}
double Erlang::sd(void){
	return sqrt(var());
}
double Erlang::f(int x){
	return f((double)x);
}
double Erlang::f(double x){
	if(x<=0.0f) return 0.0f;
	int i;
	long double sum=1.0f;
	for(i=0; i<_b; i++) sum*=_a;
	for(i=0; i<_b-1; i++) sum*=x;
	return sum*exp(-1.0f*_a*x)/_gamma;
}
double Erlang::F(int x){
	return F((double)x);
}
double Erlang::F(double x){
	if(x<=0.0f) return 0.0f;
	int i;
	double F=1.0f-1.0f/exp(_a*x), sum=1.0f/exp(_a*x);
	for(i=1; i<_b; i++){
		sum=sum*_a*x/i;
		F-=sum;
	}
	return F;
}
double Erlang::x(void){
	Exponential e(1.0f/_a);
	int i;
	double sum=0.0f;
	for(i=0; i<_b; i++) sum+=e.x();
	return sum;
}

// Chi-square Distribution

Chi_Square::Chi_Square(void): DISTRIBUTION(1.0f){
	_gamma=gamma(0.5f);
}
Chi_Square::Chi_Square(double freedom): DISTRIBUTION((int)freedom){
	if(freedom<1.0f) {
		cerr<< "warning: parameter of Chi-square distribution <1.\n";
		if(_a<0.0f) _a=-_a;
		if(_a<1.0f) _a=1.0f;
	}
	_gamma=gamma(_a/2);
}
void Chi_Square::set_a(double freedom){
	_a=(int)freedom;
	if(_a<=0) _a=-_a;
	_gamma=gamma(_a/2);
}
double Chi_Square::mean(void){
	return _a;
}
double Chi_Square::var(void){
	return 2.0f*_a;
}
double Chi_Square::sd(void){
	return sqrt(var());
}
double Chi_Square::f(int x){
	return f((double)x);
}
double Chi_Square::f(double x){
	unsigned long double output=1.0f/_gamma;
	return output*exp(-1.0f*x/2.0f)/pow(2, _a/2)*pow(x, _a/2-1);
}
double Chi_Square::F(int x){
	return F((double)x);
}
double Chi_Square::F(double x){
	if(x<=0) return 0.0f;
	else{
		double delta=x/5000.0f;
		double output=4.0f*delta*f(delta)/3.0f;
		double k, d, ave=mean();
		for(k=2.0f*delta; k+2.0f*delta<x; k+=2.0f*delta){
			output+=2.0f*delta*f(k)/3.0f;
			d=4.0f*delta*f(k+delta)/3.0f;
			output+=d;
			if(d<ERR && k>ave) break;
		}
		return output+(x-k)*f(x)/3.0f;
	}
}
double Chi_Square::x(void){
	Normal N(0, 1);
	double S2;	//sample variance
	double *z=new double [ ( unsigned int )_a + 1];
	int i;
	for(i=0; i<_a+1; i++) z[i]=N.x();
	S2=Var(z, ( unsigned int )_a + 1 );
	delete z;
	return S2*_a;
}

// T-Distribution
// _gamma= gamma[(k+1)/2]/ gamma(k/2)

T_Distribution::T_Distribution(void): DISTRIBUTION(3.0f){
	_gamma=0.664670194f;	// gamma(2.5)/gamma(2.0)
}
T_Distribution::T_Distribution(double freedom): DISTRIBUTION((int)freedom){
	if(freedom<1.0f) {
		cerr<< "warning: parameter of t-distribution <1.\n";
		if(_a<0.0f) _a=-_a;
		if(_a<1.0f) _a=1.0f;
	}
	_gamma=gamma((_a+1.0f)/2.0f)/gamma(_a/2.0f);
}
void T_Distribution::set_a(double freedom){
	_a=(int)freedom;
	if(_a<0) _a=-_a;
	if(_a<1.0f) _a=1.0f;
	_gamma=gamma((_a+1.0f)/2.0f)/gamma(_a/2.0f);
}
double T_Distribution::mean(void){
	return 0.0f;
}
double T_Distribution::var(void){
	if(_a>2.0f) return _a/(_a-2.0f);
	else if (_a>1.0f) return tan(M_PI/2.0f);
	else {
		cerr<< "warning: t-distribution, variance dosn't defined.\n";
		return 0.0f;
	}
}
double T_Distribution::sd(void){
	if(_a>2.0f) return sqrt(var());
	else if (_a>1.0f) return tan(M_PI/2.0f);
	else {
		cerr<< "warning: t-distribution, standart deviation dosn't defined.\n";
		return 0.0f;
	}
}
double T_Distribution::f(int x){
	return f((double)x);
}
double T_Distribution::f(double x){
	unsigned long double output=_gamma/sqrt(M_PI*_a);
	return output/pow(x*x/_a+1.0f, (_a+1.0f)/2.0f);
}
double T_Distribution::F(int x){
	return F((double)x);
}
double T_Distribution::F(double x){
	if(x==0.0f) return 0.5f;
	double z=fabs(x);
	double delta=z/5000.0f;
	double output=(delta*f(0)+4.0f*delta*f(delta))/3.0f;
	double k, d, ave=mean();
	for(k=2.0f*delta; k+2.0f*delta<z; k+=2.0f*delta){
		output+=2.0f*delta*f(k)/3.0f;
		d=4.0f*delta*f(k+delta)/3.0f;
		output+=d;
		if(d<ERR) break;
	}
	output+=(z-k)*f(x)/3.0f;
	if(x>0) return output+0.5f;
	else return 0.5f-output;
}
double T_Distribution::x(void){
	Normal N(0, 1);
	double X, S2;	//sample mean and variance
	double *z=new double [ ( unsigned int )_a+1];
	int i;
	for(i=0; i<_a+1; i++) z[i]=N.x();
	X=Ave(z, ( unsigned int )_a+1);
	S2=Var(z, ( unsigned int )_a+1);
	delete z;
	return X/sqrt(S2/(_a+1.0f));
}

// F-Distribution
// _gamma= gamma[(k+1)/2]/ _gamma2: gamma(k/2)

F_Distribution::F_Distribution(void): DISTRIBUTION(1, 1){
	_beta=0.101321183f;	
}
F_Distribution::F_Distribution(double u, double v): DISTRIBUTION((int)u, (int)v){
	if(u<1.0f || v<1.0f) {
		cerr<< "warning: parameter of f-distribution <1.\n";
		if(_a<0.0f) _a=-_a;
		if(_b<0.0f) _b=-_b;
		if(_a<1.0f) _a=1.0f;
		if(_b<1.0f) _b=1.0f;
	}
	_beta=Beta(_a/2.0f, _b/2.0f);
}
void F_Distribution::set_a(double u){
	_a=(int)u;
	if(_a<0) _a=-_a;
	if(_a<1.0f) _a=1.0f;
	_beta=Beta(_a/2.0f, _b/2.0f);
}
void F_Distribution::set_b(double v){
	_b=(int)v;
	if(_b<0) _b=-_b;
	if(_b<1.0f) _b=1.0f;
	_beta=Beta(_a/2.0f, _b/2.0f);
}
double F_Distribution::mean(void){
	if(_b>2.0f) return _b/(_b-2.0f);
	else {
		cerr<< "warning: f-distribution, mean dosn't defined.\n";
		return 0.0f;
	}
}
double F_Distribution::var(void){
	if(_b>4.0f) return 2.0f*_b*_b*(_a+_b-2)/(_a*(_b-2)*(_b-2)*(_b-4));
	else {
		cerr<< "warning: f-distribution, variance dosn't defined.\n";
		return 0.0f;
	}
}
double F_Distribution::sd(void){
	if(_b>4.0f) return sqrt(var());
	else {
		cerr<< "warning: f-distribution, standart deviation dosn't defined.\n";
		return 0.0f;
	}
}
double F_Distribution::f(int x){
	return f((double)x);
}
double F_Distribution::f(double x){
	unsigned long double output=1.0f/_beta;
	output*=pow(_a/_b, _a/2.0f);
	output*=pow(x, (_a/2.0f)-1.0f);
	return output/pow(x*_a/_b+1.0f, (_a+_b)/2.0f);
}
double F_Distribution::F(int x){
	return F((double)x);
}
double F_Distribution::F(double x){
	if(x<=0) return 0.0f;
	else{
		double delta=x/5000.0f;
		double output=4.0f*delta*f(delta)/3.0f;
		double k, d, ave=mean();
		for(k=2.0f*delta; k+2.0f*delta<x; k+=2.0f*delta){
			output+=2.0f*delta*f(k)/3.0f;
			d=4.0f*delta*f(k+delta)/3.0f;
			output+=d;
			if(d<ERR && k>ave) break;
		}
		return output+(x-k)*f(x)/3.0f;
	}
}
double F_Distribution::x(void){
	Normal N(0, 1);
	int i;
	double X, S2;	//sample mean and variance
	double *z=new double [ (unsigned int )_a+1];
	for(i=0; i<_a+1; i++) z[i]=N.x();
	X=Ave(z, ( unsigned int )_a+1);
	S2=Var(z, ( unsigned int )_a+1);
	delete z;
	return X/sqrt(S2/(_a+1.0f));
}
