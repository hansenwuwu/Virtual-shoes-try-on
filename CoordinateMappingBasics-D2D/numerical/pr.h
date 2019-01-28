#ifndef _PATTERN_RECOGNITION_H_
#define _PATTERN_RECOGNITION_H_

#include <numerical\distribution.h>
#include <stdarg.h>

//
// typedef ans structrue
//

typedef float		Probf, Classf;
typedef double		Probd, Classd;
typedef float *		Featuref;
typedef double *	Featured;
typedef Featuref *	FeatureArrayf;
typedef Featured *	FeatureArrayd;
typedef Classf *	ClassArrayf;
typedef Classd *	ClassArrayd;
typedef int			FeatureDim, ArrayNum;
typedef float ( *LearningRatef )( int );
typedef double ( *LearningRated )( int );

typedef struct _feature_set_float{
	FeatureArrayf x;
	ArrayNum n;
} FeatureSetf, TrueFeatureSetf, FalseFeatureSetf;

typedef struct _linear_classifier_float{
	float *w;
	FeatureDim d;
} LClassifierf;

typedef struct _linear_classifier_double{
	double *w;
	FeatureDim d;
} LClassifier, LClassifierd;

typedef struct _2lp_model{
	double *w;
	FeatureDim d;
	int k;
} NeuralNet;

typedef struct _adaboost{
	double *w, *a;
	FeatureDim d;
	int k;
} AdaBoost;

#define NewFeature( d, type )		( type * )malloc( d * sizeof( type ) )
#define NewFeatureArray( n, type )	( type ** )malloc( n * sizeof( type* ) )
#define NewClassArray( n, type )	( type * )malloc( n * sizeof( type ) )

//
// Structure handling
//

void NewLClassifierf( LClassifierf *, FeatureDim d );
void DeleteLClasifierf( LClassifierf );

void NewLClassifier( LClassifier *, FeatureDim d );
void NewNeuralNet( NeuralNet *, FeatureDim d, int k );
void NewAdaBoost( AdaBoost *, FeatureDim d, int k );
void DeleteLClassifier( LClassifier );
void DeleteNeuralNet( NeuralNet );
void DeleteAdaBoost( AdaBoost );

//
// Classification
//

int Bayesian( distribution::Normal_ND [], double [], int c );
int Euclidean( distribution::Normal_ND [], double [], int c );
int Mahalanobis( distribution::Normal_ND [], double [], int c );

// Linear Classifier

void Perceptronf( LClassifierf *, FeatureArrayf, ClassArrayf, ArrayNum, LearningRatef, int time );
void Perceptronf( LClassifierf *, TrueFeatureSetf, FalseFeatureSetf, LearningRatef, int time );
void LMS_Classifyf( LClassifierf *, FeatureArrayf, ClassArrayf, ArrayNum, LearningRatef, int time );
void LMS_Classifyf( LClassifierf *, TrueFeatureSetf, FalseFeatureSetf, LearningRatef, int time );
void LS_Classifyf( LClassifierf *, FeatureArrayf, ClassArrayf, ArrayNum );
void LS_Classifyf( LClassifierf *, TrueFeatureSetf, FalseFeatureSetf );
float Classifyf( Featuref, LClassifierf );

void Perceptron( double x[], double y[], LClassifier *, int nx, double lr, int time = 50 );
void LMS_Train( double x[], double y[], LClassifier *, int nx, double lr );
void LS_Train( double x[], double y[], LClassifier *, int nx );
double Classify( double x[], LClassifier );

// Nonlinear Classifier

double BackProp( double x[], double y[], NeuralNet *, int n, 
			  double lr, double ri, double rd, double c, int time );
void Boosting( double x[], double y[], AdaBoost *, int n );
double Classify( double x[], NeuralNet );
double Classify( double x[], AdaBoost );

//
// Estimation
//

// Expectation-Maximization
void EM_Estimate( double x[], int nx, int d, int n, distribution::Mixture *, int time = 50 );

// Parsen Window
double ParzenWin( double xi[], double x[], int nx, int d, double h );

//
// Feature Selection
//

void ScatterB( double *sb, double *x, int m, int n, int d );
void ScatterW( double *sw, double *x, int m, int n, int d );
void ScatterM( double *sm, double *x, int m, int n, int d );
void NPScatterB( double *sb, double *x, int m, int n, int d, int knn, double a );
double J1( double *x, int m, int n, int d );
double J2( double *x, int m, int n, int d );
double J3( double *x, int m, int n, int d );
double Fdr( double *x, int m, int n, int d );
/*
	s: d*d matrix
	xmn is a d*1 vector
	m class
	n data within a class
	x = x11, ..., x1n, x21, ..., x2n, ... xm1, ..., xmn
*/

double J1( int d, int m, ... );
double J2( int d, int m, ... );
double J3( int d, int m, ... );
double Fdr( int d, int m, ... );
/*	
	example: J3( d, m, x1, ..., xm, n1, ..., nm );
	x1, ..., xm: data
	n1, ..., nm: number of data
*/

//
// Feature Generation
//

int pca( double *pc, double *var, double *x, int n, int d );
int pca( double *pc, double *x, int n, int d, int *l, const double r );
/*	pc: d * d matrix A -> d * l matrix A	*/
int pca( double *pc, double *x, int n, int d, int l );
/*	pc: d * l matrix A, y = ATx				*/

int pda( double *a, double *x, int m, int n, int d );
int nda( double *a, double *x, int m, int n, int d, int k, double param );
/*	a: d * d matrix, y = aTx				*/

#endif