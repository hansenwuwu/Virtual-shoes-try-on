#ifndef _STATISTIC_H_
#define _STATISTIC_H_

#define SQUARE(x) (x)*(x)

// Statistic functions

double Ave(int x[], int n);
double Ave(double x[], int n);
double Var(int x[], int n);
double Var(double x[], int n);
double Cov(int x[], int y[], int n);
double Cov(double x[], double y[], int n);
double Corr(int x[], int y[], int n);
double Corr(double x[], double y[], int n);

// Regration function

double Sxx(int x[], int n);
double Sxx(double x[], int n);
double Sxy(int x[], int y[], int n);
double Sxy(double x[], double y[], int n);
double beta1(int x[], int y[], int n);
double beta1(double x[], double y[], int n);
double beta0(int x[], int y[], int n);
double beta0(double x[], double y[], int n);
double SSE(int x[], int y[], int n);
double SSE(double x[], double y[], int n);
double SSR(int x[], int y[], int n);
double SSR(double x[], double y[], int n);
double R2(int x[], int y[], int n);
double R2(double x[], double y[], int n);

// Confidence Interval

#endif