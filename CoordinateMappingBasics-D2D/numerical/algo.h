#ifndef _ALGO_H_
#define _ALGO_H_

#include <iostream>
using namespace std;

// Clustering

/*
K-mean Clustering
	data[n * d]
	cluster[c * d]
	domain[2 * d] = [ min0, max0, min1, max1, ...]
	out[n]
*/
void Kmean( int *data, int *cluster, int *domain, int *output, int n, int d, int c, 
		   bool init, void *mem = 0x00 ); 
void Kmean( double *data, double *cluster, double *domain, int *output, int n, int d, int c, 
		   bool init, void *mem = 0x00 ); 

// Histogram

void Histogram( int *data, int bottom, int top, int n, int *H = NULL );
void Histogram( int *source, int *target, int b, int t, 
			   int n1, int n2, int *H1 = NULL, int *H2 = NULL );
void Histogram( int *source, int *target, int *bottom, int *top, int d, 
			   int n1, int n2, int time );

// Normalized Cut

double NormalizedCut( double *W, double *x, int n, int time );

// Sorting

template < typename T > void BubbleSort( T *a, int n );
template < typename T > void SelectSort( T *a, int n );
template < typename T > void InsertSort( T *sort, int n );
template < typename T > void QuickSort( T *a, int n );

//
// Implement
//

template < typename T > void BubbleSort( T *a, int n ){
	int i, j;
	for( i = 0; i < n - 1; i++ ){
		for( j = 0; j < n - i; j++ )
			if( a[j] < a[j + 1] ) swap( a[j], a[j + 1] );
	}
}

template < typename T > void SelectSort( T *a, int n ){
	int i, j, flag;
	for( i = 0; i < n - 1; i++ ){
		flag = i;
		for( j = i + 1; j < n; j++ )
			if( a[j] < a[flag] ) flag = j;
		swap( flag, j );
	}
}

template < typename T > void InsertSort( T *a, int n ){
	if( n == 1 ) return;

	T temp;
	int i, j, k;
	
	if( a[0] > a[1] ) swap( a[0], a[1] );
	for( i = 2; i < n; i++ ){			// a[i] need to be insert
		temp = a[i];
		for( j = i - 1; j >= 0; j-- ){
			if( a[j] <= a[i] ){
				for( k = i; k - 1 > j; k-- ) a[k] = a[k - 1];
				a[k] = temp;
				break;
			}
		}
		if( j < 0 ){				// smaller then all
			for( k = i; k - 1 >= 0; k-- ) a[k] = a[k - 1];
			a[k] = temp;
		}
	}
}

template < typename T > void QuickSort( T *a, int n ){
	if( n < 1 ) return;
	int i = 1, j = n - 1;

	swap( a[0], a[n / 2] );
	if( i < j ){
		while( i < j ){
			while( a[i] <= a[0] && i < n ) i++;
			while( a[j] >= a[0] && j > 0 ) j--;
			if( i < j ) swap( a[i], a[j] );
		}
		swap( a[0], a[j] );
	}

	QuickSort( a, j );
	QuickSort( a + j + 1, n - j - 1 );
}

#endif