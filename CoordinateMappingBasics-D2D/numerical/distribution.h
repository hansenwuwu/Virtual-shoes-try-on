#ifndef _DISTRIBUTION_H_
#define _DISTRIBUTION_H_

namespace distribution{

//
// 1D distribution class
//

class DISTRIBUTION{
protected:
	double _a, _b;
public:
	DISTRIBUTION(void);
	DISTRIBUTION(double a);
	DISTRIBUTION(double a, double b);

	virtual void set_a(double a);
	virtual void set_b(double b);
	virtual double mean(void)=0;	//E(X)
	virtual double var(void)=0;		//V(X)
	virtual double sd(void)=0;		//s.d., standart deviation
	virtual double f(int x)=0;		//f(x) for discrete distribution
	virtual double f(double x)=0;	//f(x) for contineous distribution
	virtual double F(int x)=0;		//F(x) for discrete distribution
	virtual double F(double x)=0;	//F(x) for contineous distribution
	virtual double x(void)=0;		//generate a random value follow this distribution
};
typedef DISTRIBUTION Distribution;

//
// Discrete Random Variables
//

// Discrete Uniform distribution: 
// f(x)=1/n;
class Uniform_Disc: public DISTRIBUTION{
public:
	Uniform_Disc(void);
	Uniform_Disc(double a, double b);

	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// Poisson distribution: 
// f(x)=(e^(-a)*a^x)/x!
class Poisson: public DISTRIBUTION{
public:
	Poisson(void);
	Poisson(double ave);

	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

//
// Continuous Random Variables
//

// Continuous Uniform distribution: 
// f(x)=1/n;
class Uniform_Cont: public DISTRIBUTION{
public:
	Uniform_Cont(void);
	Uniform_Cont(double a, double b);

	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// Normal distribution: 
// f(x)=1/(sqrt(2pi)*b)*e^(-0.5*(x-a)^2/b^2);
class Normal: public DISTRIBUTION{
public:
	Normal(void);	//standard normal distribution
	Normal(double ave, double var);

	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// Exponential distribution:
// f(x)=a*e^(-ax)
class Exponential: public DISTRIBUTION{
public:
	Exponential(void);
	Exponential(double ave);	//input: average of distribution

	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// Gamma distribution:
// f(x)=(a^b*x^(b-1)*e^-ax)/gamma(b);
// a=lambda, b=r.
class Gamma: public DISTRIBUTION{
protected:
	unsigned long double _gamma;
public:
	Gamma(void);
	Gamma(double scale, double shape);	//input: scale and shape parameters

	virtual void set_b(double shape);	//input: the shape parameter
	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// Erlang distribution:
// f(x)=(a^b*x^(b-1)*e^-ax)/gamma(b);
// a=lambda, b=r=1, 2, 3, ....
class Erlang: public DISTRIBUTION{
protected:
	unsigned long double _gamma;
public:
	Erlang(void);
	Erlang(double scale, double shape);	//input: scale and shape parameters

	virtual void set_b(double shape);	//input: the shape parameter
	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// Chi-square distribution:
// f(x)=1/(2^(a/2)*gamma(k/2))*x^(a/2-1)*e^(-x/2);
// a: degree of freedom, must be a positive interger
class Chi_Square: public DISTRIBUTION{
protected:
	unsigned long double _gamma;
public:
	Chi_Square(void);
	Chi_Square(double freedom);	//input: degrees of freedom

	virtual void set_a(double freedom);	//input: degrees of freedom
	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// T distribution:
// f(x)=[gamma((a+1)/2)/(sqrt(a*pi)*gamma(a/2))]/[(x^2/a)+1]^((a+1)/2)
// a: degree of freedom, must be a positive interger
class T_Distribution: public DISTRIBUTION{
protected:
	unsigned long double _gamma;
public:
	T_Distribution(void);
	T_Distribution(double freedom);	//input: degrees of freedom

	virtual void set_a(double freedom);	//input: degrees of freedom
	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

// F-distribution:
// f(x)=(a/b)^(a/2)*x^((a/2)-1)/[beta(a/2, b/2)*((a/b)x+1)^((a+b)/2)];
// a=lambda, b=r.
class F_Distribution: public DISTRIBUTION{
protected:
	unsigned long double _beta;
public:
	F_Distribution(void);
	F_Distribution(double u, double v);	
	//input: degrees of freedom in the numerator and denominator

	virtual void set_a(double u);
	virtual void set_b(double v);	//input: the shape parameter
	virtual double mean(void);	//E(X)
	virtual double var(void);		//V(X)
	virtual double sd(void);		//s.d., standart deviation
	virtual double f(int x);		//f(x) for discrete distribution
	virtual double f(double x);	//f(x) for contineous distribution
	virtual double F(int x);		//F(x) for discrete distribution
	virtual double F(double x);	//F(x) for contineous distribution
	virtual double x(void);		//generate a random value follow this distribution
};

//
// N-D Distribution Basic Class
//

class Distribution_ND{
protected:
	int _d, _m_size;
	double *_a, *_b;
	void *_mem;
public:
	Distribution_ND( void );
	Distribution_ND( int d );
	Distribution_ND( int d, double *a, double *b );
	~Distribution_ND( void );

	int GetDimension( void );
	virtual void SetMean( double * ) = 0;
	virtual void SetVar( double * ) = 0;
	virtual void Mean( double * ) = 0;	// E(X)
	virtual void Var( double * ) = 0;	// V(X)	
	virtual double f( double * ) = 0;
	virtual void Rand( double * ) = 0;		// generate a random value follow this distribution
};

class Normal_ND: public Distribution_ND{
public:
	Normal_ND( int d );
	Normal_ND( int d, double *a, double *b );

	virtual void SetMean( double * );
	virtual void SetVar( double * );
	virtual void Mean( double * );		// E(X)
	virtual void Var( double * );		// V(X)
	virtual double f( double * );
	virtual void Rand( double * );			// generate a random value follow this distribution
};

class Mixture: public Distribution_ND{
protected:
	Distribution_ND **_X;

public:
	Mixture( int n, Distribution_ND * );
	~Mixture( void );

	void SetProb( double * );
	void SetMean( int, double * );
	void SetVar( int, double * );
	double Prob( int );
	void Mean( int, double * );
	void Var( int, double * );
	double f( int, double * );

	virtual void SetMean( double * );
	virtual void SetVar( double * );
	virtual void Mean( double * );
	virtual void Var( double * );
	virtual double f( double * );
	virtual void Rand( double * );
};

}	// namespace distribution;

#endif