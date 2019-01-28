#include <math.h>
#include "Statistic.h"
#include "distribution.h"

// inline function

template< typename T > T SQR( T x ){
	return x * x;
}

//
// Statistic functions
//

double Ave(int* x, int n){
	long int sum;
	int i=0;
	for(i=0, sum=0; i<n; i++) sum+=x[i];
	return (double)sum/n;
};
double Ave(double* x, int n){
	double sum=0;
	for(int i=0; i<n; i++) sum+=x[i];
	return sum/n;
};

double Var(int* x, int n){
	double ave=Ave(x, n), sum=0;
	for(unsigned long i=0; i<(unsigned)n; i++) sum+=(x[i]-ave)*(x[i]-ave);
	return sum/(n-1);
};
double Var(double* x, int n){
	double ave=Ave(x, n);
	unsigned long double sum=0;
	for(unsigned long i=0; i<(unsigned)n; i++) sum+=(x[i]-ave)*(x[i]-ave);
	return sum/(n-1);
};

double Cov(int* x, int* y, int n){
	int i;
	long double sum, x_ave=Ave(x, n);
	for(i=0, sum=0; i<n; i++) sum+=y[i]*(x[i]-x_ave);
	return (double)sum/(n-1);
}
double Cov(double* x, double* y, int n){
	int i;
	long double sum, x_ave=Ave(x, n);
	for(i=0, sum=0; i<n; i++) sum+=y[i]*(x[i]-x_ave);
	return (double)sum/(n-1);
}

double Corr(int* x, int* y, int n){
	return Cov(x, y, n)/sqrt(Var(x, n)*Var(y, n));
}
double Corr(double* x, double* y, int n){
	return Cov(x, y, n)/sqrt(Var(x, n)*Var(y, n));
}

//
// Regration function
//

double Sxx(int* x, int n){
	unsigned long double sum=0; 
	double ave=Ave(x, n);
	for(int i=0; i<n; i++) sum+=(x[i]-ave)*(x[i]-ave);
	return sum;
}
double Sxx(double* x, int n){
	unsigned long double sum=0, ave=Ave(x, n);
	for(int i=0; i<n; i++) sum+=(x[i]-ave)*(x[i]-ave);
	return sum;
}
double Sxy(int* x, int* y, int n){
	long double sum=0;
	double x_ave=Ave(x, n), y_ave=Ave(x, n);
	for(int i=0; i<n; i++) sum+=(x[i]-x_ave)*(y[i]-y_ave);
	return sum;
}
double Sxy(double* x, double* y, int n){
	long double sum=0;
	double x_ave=Ave(x, n), y_ave=Ave(x, n);
	for(int i=0; i<n; i++) sum+=(x[i]-x_ave)*(y[i]-y_ave);
	return sum;
}
double beta1(int* x, int* y, int n){
	return Sxy(x, y, n)/Sxx(x, n);
}
double beta1(double* x, double* y, int n){
	return Sxy(x, y, n)/Sxx(x, n);
}
double beta0(int* x, int* y, int n){
	return Ave(y, n)-beta1(x, y, n)*Ave(x, n);
}
double beta0(double* x, double* y, int n){
	return Ave(y, n)-beta1(x, y, n)*Ave(x, n);
}
double SSE(int* x, int* y, int n){
	return Sxx(y, n)-beta1(x, y, n)*Sxy(x, y, n);
}
double SSE(double* x, double* y, int n){
	return Sxx(y, n)-beta1(x, y, n)*Sxy(x, y, n);
}
double SSR(int* x, int* y, int n){
	return Sxx(y, n)-SSE(x, y, n);
}
double SSR(double* x, double* y, int n){
	return Sxx(y, n)-SSE(x, y, n);
}
double R2(int* x, int* y, int n){
	return SSR(x, y, n)/Sxx(y, n);
}
double R2(double* x, double* y, int n){
	return SSR(x, y, n)/Sxx(y, n);
}
