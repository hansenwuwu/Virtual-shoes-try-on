 
#include "track.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

// k-dtree
#include "kdtree.h"

//
// thread global variables
//

__declspec( thread )  extern double gxd, gyd, gzd, gcad, gcbd, gccd, gsad, gsbd, gscd;
__declspec( thread )  extern float gxf, gyf, gzf, gcaf, gcbf, gccf, gsaf, gsbf, gscf;

// 
// macro
//

#define M_PI			3.1415926535897932384626433832795
#define SQR(x)			( (x) * (x) )
#define CompDist2f( d2, m, p1, p2 ) \
	gxf = (m)[0][0] * (p1).x + (m)[0][1] * (p1).y + (m)[0][2] * (p1).z + (m)[0][3] - (p2).x;\
	gyf = (m)[1][0] * (p1).x + (m)[1][1] * (p1).y + (m)[1][2] * (p1).z + (m)[1][3] - (p2).y;\
	gzf = (m)[2][0] * (p1).x + (m)[2][1] * (p1).y + (m)[2][2] * (p1).z + (m)[2][3] - (p2).z;\
	(d2) = SQR( gxf ) + SQR( gyf ) + SQR( gzf );

#define CompDist2d( d2, m, p1, p2 ) \
	gxd = (m)[0][0] * (p1).x + (m)[0][1] * (p1).y + (m)[0][2] * (p1).z + (m)[0][3] - (p2).x;\
	gyd = (m)[1][0] * (p1).x + (m)[1][1] * (p1).y + (m)[1][2] * (p1).z + (m)[1][3] - (p2).y;\
	gzd = (m)[2][0] * (p1).x + (m)[2][1] * (p1).y + (m)[2][2] * (p1).z + (m)[2][3] - (p2).z;\
	(d2) = SQR( gxd ) + SQR( gyd ) + SQR( gzd );

#define CompDist2fv( d2, m, p1, p2 ) \
	gxf = (m)[0][0] * (p1).ptr[0] + (m)[0][1] * (p1).ptr[1] + (m)[0][2] * (p1).ptr[2] + (m)[0][3] - (p2).ptr[0];\
	gyf = (m)[1][0] * (p1).ptr[0] + (m)[1][1] * (p1).ptr[1] + (m)[1][2] * (p1).ptr[2] + (m)[1][3] - (p2).ptr[1];\
	gzf = (m)[2][0] * (p1).ptr[0] + (m)[2][1] * (p1).ptr[1] + (m)[2][2] * (p1).ptr[2] + (m)[2][3] - (p2).ptr[2];\
	(d2) = SQR( gxf ) + SQR( gyf ) + SQR( gzf );

#define CompDist2dv( d2, m, p1, p2 ) \
	gxd = (m)[0][0] * (p1).ptr[0] + (m)[0][1] * (p1).ptr[1] + (m)[0][2] * (p1).ptr[2] + (m)[0][3] - (p2).ptr[0];\
	gyd = (m)[1][0] * (p1).ptr[0] + (m)[1][1] * (p1).ptr[1] + (m)[1][2] * (p1).ptr[2] + (m)[1][3] - (p2).ptr[1];\
	gzd = (m)[2][0] * (p1).ptr[0] + (m)[2][1] * (p1).ptr[1] + (m)[2][2] * (p1).ptr[2] + (m)[2][3] - (p2).ptr[2];\
	(d2) = SQR( gxd ) + SQR( gyd ) + SQR( gzd );

#define CompDistvf( v, m, p1, p2 ) \
	(v).x = (m)[0][0] * (p1).x + (m)[0][1] * (p1).y + (m)[0][2] * (p1).z + (m)[0][3] - (p2).x;\
	(v).y= (m)[1][0] * (p1).x + (m)[1][1] * (p1).y + (m)[1][2] * (p1).z + (m)[1][3] - (p2).y;\
	(v).z = (m)[2][0] * (p1).x + (m)[2][1] * (p1).y + (m)[2][2] * (p1).z + (m)[2][3] - (p2).z;

#define NewPoints( v, m, p1 ) \
	(v).x[0] = (m)[0][0] * (p1).x + (m)[0][1] * (p1).y + (m)[0][2] * (p1).z + (m)[0][3]; \
	(v).x[1] = (m)[1][0] * (p1).x + (m)[1][1] * (p1).y + (m)[1][2] * (p1).z + (m)[1][3]; \
	(v).x[2] = (m)[2][0] * (p1).x + (m)[2][1] * (p1).y + (m)[2][2] * (p1).z + (m)[2][3];

#define CompDistvf_kd( v, m, p1, p2 ) \
	(v).x = (m)[0][0] * (p1).x + (m)[0][1] * (p1).y + (m)[0][2] * (p1).z + (m)[0][3] - (p2).x; \
	(v).y = (m)[1][0] * (p1).x + (m)[1][1] * (p1).y + (m)[1][2] * (p1).z + (m)[1][3] - (p2).y; \
	(v).z = (m)[2][0] * (p1).x + (m)[2][1] * (p1).y + (m)[2][2] * (p1).z + (m)[2][3] - (p2).z;


void swap_float(float *one, float *two)
{
	float temp_swap;
	temp_swap = *one;
	*one = *two;
	*two = temp_swap;
}

void swap_Point3f(Point3f *one, Point3f *two)
{
	Point3f temp_swap;
	temp_swap.x = one->x;
	temp_swap.y = one->y;
	temp_swap.z = one->z;
	one->x = two->x;
	one->y = two->y;
	one->z = two->z;
	two->x = temp_swap.x;
	two->y = temp_swap.y;
	two->z = temp_swap.z;
}

void QuickSortICP(Point3f buf[],float d[], int l, int h)
{
	int i,j;
	float key;
	Point3f key_buf;
	if (l>=h)return ;
	i = l; j = h; 
	key = d[i];
	key_buf.x = buf[i].x;
	key_buf.y = buf[i].y;
	key_buf.z = buf[i].z;
	while( i<j )
	{
		while( i<j && d[j]>key) j--;
		if (i<j)
		{
			d[i] = d[j];
			buf[i].x = buf[j].x;
			buf[i].y = buf[j].y;
			buf[i].z = buf[j].z;
			i++;
		}
		while ( i<j && d[i]<key )i++;
		if (i<j)
		{
			d[j] = d[i];
			buf[j].x = buf[i].x;
			buf[j].y = buf[i].y;
			buf[j].z = buf[i].z;
			j--;
		}
	}
	d[i] = key;
	buf[i].x = key_buf.x;
	buf[i].y = key_buf.y;
	buf[i].z = key_buf.z;
	if (l<i-1)
		QuickSortICP(buf,d,l,i-1);
	if (i+1<h)
		QuickSortICP(buf,d,i+1,h);

}

//
// Data Transform by Points by ICP
//

float ICPf( Point3f_Array patch, Point3f_Array base, RT_Paramf *theta, const ICP_Paramf param, const Point3f_Array buf , float icpe[100] ){
	float m[3][4], e, d2;
	struct _point_3d_float *pp, *pbase, *pbuf, *np;
	struct _point_3d_float_pair_array pair;
	Point3f temp, vec;
	int i, j, k;

	// initialize
	pair.pt1 = patch.pt;	// point of patch
	pair.pt2 = buf.pt;		// nearest point on base
	pair.n = patch.n;

	for( k = 0; k < param.iter; k++ ){
		memset( &vec, 0x00, sizeof( Point3f ) );
		// update transform matrix
		CompRTMatf( m, *theta );

		// update points's pair
		pp = patch.pt;
		pbuf = buf.pt;
		for( i = 0; i < patch.n; i++, pp++, pbuf++ ){	// each point on patch 
			pbase = base.pt;
			e = FLT_MAX;
			// old method
			for( j = 0; j < base.n; j++, pbase++ ){		// find nearset point on base
				CompDistvf( temp, m, *pp, *pbase );		// compute distance from patch i to base j
				d2 = SQR( temp.x ) + SQR( temp.y ) + SQR( temp.z );
				if( d2 < e ){							// if is minimal
					*pbuf = *pbase;						// save data
					e = d2;							// update minimal distance
				}
			}

			CompDistvf( temp, m, *pp, *pbuf );
			vec.x += temp.x;	vec.y += temp.y;	vec.z += temp.z;
		}
		// Rigid body transform by SVD
		RigidTransform_SVDMatf( pair, m );
		CompRTFromMatf( theta, m );
		// translation
		theta->x -= vec.x / ( float )patch.n;
		theta->y -= vec.y / ( float )patch.n;
		theta->z -= vec.z / ( float )patch.n;
		// transformation
		e = RigidTransformf( pair, theta, param.trans );
		e = e/patch.n;
		//e = 0.5;
		// normalize the rotate angle
		theta->a = ( float )fmod( ( float )theta->a, ( float )M_PI );
		theta->b = ( float )fmod( ( float )theta->b, ( float )M_PI );
		theta->c = ( float )fmod( ( float )theta->c, ( float )M_PI );
		// record error
		icpe[k+1] = e;
		icpe[0] = k+1;
		if( e < param.trans.err ) break;
	}
	return e;
}

float ICP_kdtree(Point3f_Array patch, Point3f_Array base, RT_Paramf *theta, const ICP_Paramf param, const Point3f_Array buf, float icpe[100]){
	float m[3][4], e, d2;
	struct _point_3d_float *pp, *pbase, *pbuf, *np, *pbuf1;
	struct _point_3d_float_pair_array pair;
	Point3f temp, vec;
	int i, j, k;
	// k-d tree
	struct kd_node_t *kd, *root, *found, data;
	double best_dist;
	// time
	unsigned int start, end;

	// construct k-d tree
	pbase = base.pt;
	int N = base.n; // The number of nodes
	kd = (struct kd_node_t*) calloc(N, sizeof(struct kd_node_t));
	for (i = 0; i < N; i++, pbase++){
		kd[i].x[0] = (double)(*pbase).x;
		kd[i].x[1] = (double)(*pbase).y;
		kd[i].x[2] = (double)(*pbase).z;
	}
	root = make_tree(kd, N, 0, 3);
	
	// initialize
	pair.pt1 = patch.pt;	// point of patch
	pair.pt2 = buf.pt;		// nearest point on base
	pair.n = patch.n;	

	for (k = 0; k < param.iter; k++){
		memset(&vec, 0x00, sizeof(Point3f));
		// update transform matrix
		CompRTMatf(m, *theta);

		// update points's pair
		pp = patch.pt;
		pbuf = buf.pt;
		for (i = 0; i < patch.n; i++, pp++, pbuf++){	// each point on patch 
			pbase = base.pt;
			e = FLT_MAX;

			// find nearset point using k-d tree
			NewPoints(data, m, *pp);
			visited = 0;
			found = 0;
			nearest(root, &data, 0, 3, &found, &best_dist);
			(*pbuf).x = found->x[0]; (*pbuf).y = found->x[1]; (*pbuf).z = found->x[2];

			CompDistvf(temp, m, *pp, *pbuf);
			vec.x += temp.x;	vec.y += temp.y;	vec.z += temp.z;
		}
		
		// Rigid body transform by SVD
		RigidTransform_SVDMatf(pair, m);
		CompRTFromMatf(theta, m);
		// translation
		theta->x -= vec.x / (float)patch.n;
		theta->y -= vec.y / (float)patch.n;
		theta->z -= vec.z / (float)patch.n;
		// transformation
		e = RigidTransformf(pair, theta, param.trans);
		e = e / patch.n;
		// normalize the rotate angle
		theta->a = (float)fmod((float)theta->a, (float)M_PI);
		theta->b = (float)fmod((float)theta->b, (float)M_PI);
		theta->c = (float)fmod((float)theta->c, (float)M_PI);
		// record error
		icpe[k + 1] = e;
		icpe[0] = k + 1;
		float err = e - icpe[k];
		if (err < 0) err = -err;
		//printf("Error = %f\n", err);
		if (err < param.trans.err) break;		
	}

	//printf("patch.n = %d, base.n = %d\n", patch.n, base.n); // 799 948
	//printf("ICP iteration times = %d\n", k);
	//printf("進行運算所花費的時間 = %.3f sec\n", (float)(end - start) / CLOCKS_PER_SEC);
	free(kd);

	return e;
}

float TrICPf( Point3f_Array patch, Point3f_Array base, RT_Paramf *theta, const ICP_Paramf param, const Point3f_Array buf , float icpe[100]){
	float m[3][4], e, d2, *d_array, swaptemp, last_e = 999;
	struct _point_3d_float *pp, *pbase, *pbuf;
	struct _point_3d_float_pair_array pair;
	Point3f temp, vec;
	int i, j, k;
	Point3f *buf_array;
	//Trimmed ICP para
 	float r = 1.0f;
 	int overlap = (int)(r * patch.n);
// 	float Slts = FLT_MAX;
	//
	d_array = (float *)malloc( patch.n*sizeof(float) );
	buf_array = (Point3f *)malloc(patch.n*sizeof(Point3f));

	

	// initialize
	pair.pt1 = patch.pt;	// point of patch
	pair.pt2 = buf.pt;		// nearest point on base
	pair.n = patch.n;

	for( k = 0; k < param.iter; k++ ){
		memset( &vec, 0x00, sizeof( Point3f ) );
		// update transform matrix
		CompRTMatf( m, *theta );

		// update points's pair
		pp = patch.pt;
		pbuf = buf.pt;
		for( i = 0; i < patch.n; i++, pp++, pbuf++ ){	// each point on patch
			pbase = base.pt;
			e = FLT_MAX;
			for( j = 0; j < base.n; j++, pbase++ ){		// find nearset point on base
				CompDistvf( temp, m, *pp, *pbase );		// compute distance from patch i to base j
				d2 = SQR( temp.x ) + SQR( temp.y ) + SQR( temp.z );
				if( d2 < e ){							// if is minimal
					*pbuf = *pbase;						// save data
					e = d2;							// update minimal distance
				}
			}
			d_array[i] = e;
			buf_array[i].x = pbuf->x;
			buf_array[i].y = pbuf->y;
			buf_array[i].z = pbuf->z;
// 			CompDistvf( temp, m, *pp, *pbuf );
// 			vec.x += temp.x;	vec.y += temp.y;	vec.z += temp.z;
		}
		
		//bubble sort
		pp = patch.pt;
		pbuf = buf.pt;

		
 		for( i = 0; i < overlap; i++, pp++, pbuf++ )
 				{
 					for ( j = i+1; j < patch.n; j++ )
 					{
 						if ( d_array[j] < d_array[i] )
 						{
 							swaptemp = d_array[j];
 							d_array[j] = d_array[i];
 							d_array[i] = swaptemp;
 							swaptemp = buf_array[j].x;
 							buf_array[j].x = buf_array[i].x;
 							buf_array[i].x = swaptemp;
 							swaptemp = buf_array[j].y;
 							buf_array[j].y = buf_array[i].y;
 							buf_array[i].y = swaptemp;
 							swaptemp = buf_array[j].z;
 							buf_array[j].z = buf_array[i].z;
 							buf_array[i].z = swaptemp;
 						}
 					}
 		
 // 			CompDistvf( temp, m, *pp,  buf_array[i] );
 // 			vec.x += temp.x;	vec.y += temp.y;	vec.z += temp.z;
 		}

		//QuickSortICP( buf_array, d_array, 0 , patch.n );
		pp = patch.pt;
		pbuf = buf.pt;
		for ( i = 0; i < overlap; i++, pp++, pbuf++ )
		{
			CompDistvf( temp, m, *pp,  buf_array[i] );
			vec.x += temp.x;	
			vec.y += temp.y;	
			vec.z += temp.z;
		}


		// Rigid body transform by SVD
		RigidTransform_SVDMatf( pair, m );
		CompRTFromMatf( theta, m );
		// translation
		theta->x -= vec.x / ( float )overlap;
		theta->y -= vec.y / ( float )overlap;
		theta->z -= vec.z / ( float )overlap;
		// transformation
		e = RigidTransformf( pair, theta, param.trans );
		e = e/( float )overlap;
		// normalize the rotate angle
		theta->a = ( float )fmod( ( float )theta->a, ( float )M_PI );
		theta->b = ( float )fmod( ( float )theta->b, ( float )M_PI );
		theta->c = ( float )fmod( ( float )theta->c, ( float )M_PI );
		
		icpe[k+1] = e;
		icpe[0] = k+1;
		if( e < param.trans.err ) break;

		if( (k > 1) && ( fabs(last_e - e ) < 0.001) )
			break;
		else
			last_e = e;
	}
	return e;
}