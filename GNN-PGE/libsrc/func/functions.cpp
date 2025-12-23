#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "func/gendef.h"




#ifdef UNIX
void strupr(char *_msg)
{
	int dist = 'A' - 'a';
	char *c_ptr = _msg;

	while (*c_ptr)
	{
		if (*c_ptr >= 'a')
			*c_ptr += dist;
		c_ptr ++;
	}
}
#endif

void error(char *t, bool ex)
{
    fprintf(stderr, t);
    if (ex)
        exit(0);
}

double area(int dimension, double *mbr)

{
    int i;
    double sum;

    sum = 1.0;
    for (i = 0; i < dimension; i++)
	sum *= mbr[2*i+1] - mbr[2*i];

    return sum;
}

double margin(int dimension, double *mbr)

{
    double *ml, *mu, *m_last, sum;

    sum = 0.0;
    m_last = mbr + 2*dimension;
    ml = mbr;
    mu = ml + 1;
    while (mu < m_last)
    {
	sum += *mu - *ml;
	ml += 2;
	mu += 2;
    }

    return sum;
}

bool inside(double &p, double &lb, double &ub)

{
    return (p >= lb && p <= ub);
}

bool inside(double *v, double *mbr, int dimension)

{
    int i;

    for (i = 0; i < dimension; i++)
	if (!inside(v[i], mbr[2*i], mbr[2*i+1]))
	    return FALSE;

    return TRUE;
}



double* overlapRect(int dimension, double *r1, double *r2)
{
        double *overlap = new double[2*dimension];
        for (int i=0; i<dimension; i++)
        {
            if ((r1[i*2]>r2[i*2+1]) || (r1[i*2+1]<r2[i*2])) 
	    {
                delete [] overlap;
		return NULL;
	    }
	    overlap[2*i] = max(r1[i*2], r2[i*2]);
            overlap[2*i+1] = min(r1[i*2+1], r2[i*2+1]);
        }

        return overlap;
}

double overlap(int dimension, double *r1, double *r2)


{
    double sum;
    double *r1pos, *r2pos, *r1last, r1_lb, r1_ub, r2_lb, r2_ub;

    sum = 1.0;
    r1pos = r1; r2pos = r2;
    r1last = r1 + 2 * dimension;

    while (r1pos < r1last)
    {
	r1_lb = *(r1pos++);
	r1_ub = *(r1pos++);
	r2_lb = *(r2pos++);
	r2_ub = *(r2pos++);

        

        if (inside(r1_ub, r2_lb, r2_ub))
        
	{
            if (inside(r1_lb, r2_lb, r2_ub))
            
                sum *= (r1_ub - r1_lb);
            else
                sum *= (r1_ub - r2_lb);
	}
	else
	{
            if (inside(r1_lb, r2_lb, r2_ub))
	    
		sum *= (r2_ub - r1_lb);
	    else 
	    {
		if (inside(r2_lb, r1_lb, r1_ub) &&
		    inside(r2_ub, r1_lb, r1_ub))
	        
		    sum *= (r2_ub - r2_lb);
		else
		
		    sum = 0.0;
	    }
	}
    }

    return sum;
}

void enlarge(int dimension, double **mbr, double *r1, double *r2)

{
    int i;

    *mbr = new double[2*dimension];
    for (i = 0; i < 2*dimension; i += 2)
    {
	(*mbr)[i]   = min(r1[i],   r2[i]);

	(*mbr)[i+1] = max(r1[i+1], r2[i+1]);
    }
}

bool section(int dimension, double *mbr1, double *mbr2)
{
    int i;

    for (i = 0; i < dimension; i++)
    {
	if (mbr1[2*i]     > mbr2[2*i + 1] ||
	    mbr1[2*i + 1] < mbr2[2*i]) 
	    return FALSE;
    }
    return TRUE;
}

int sort_lower_mbr(const void *d1, const void *d2)


{
    SortMbr *s1, *s2;
    double erg;
    int dimension;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;
    erg = s1->mbr[2*dimension] - s2->mbr[2*dimension];
    if (erg < 0.0)
	return -1;
    else if (erg == 0.0)
	return 0;
    else 
	return 1;
}

int sort_upper_mbr(const void *d1, const void *d2)


{
    SortMbr *s1, *s2;
    double erg;
    int dimension;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;
    erg = s1->mbr[2*dimension+1] - s2->mbr[2*dimension+1];
    if (erg < 0.0)
	return -1;
    else if (erg == 0.0)
	return 0;
    else 
	return 1;
}

int sort_center_mbr(const void *d1, const void *d2)

{
    SortMbr *s1, *s2;
    int i, dimension;
    double d, e1, e2;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;

    e1 = e2 = 0.0;
    for (i = 0; i < dimension; i++)
    {
        d = ((s1->mbr[2*i] + s1->mbr[2*i+1]) / 2.0) - s1->center[i];
        e1 += d*d;
        d = ((s2->mbr[2*i] + s2->mbr[2*i+1]) / 2.0) - s2->center[i];
        e2 += d*d;
    }

    if (e1 < e2)
	return -1;
    else if (e1 == e2)
	return 0;
    else 
	return 1;
}

int sortmindist(const void *element1, const void *element2)
{
    
    
    
    

	
	

    BranchList *e1,*e2;

    e1 = (BranchList *) element1;
    e2 = (BranchList *) element2;
    
    if (e1->mindist < e2->mindist) 
	return(-1);
    else if (e1->mindist > e2->mindist)
	return(1);
    else
	return(0);

}

int pruneBrunchList(double *nearest_distanz, const void *activebrunchList, 
		    int n)
{

    
    
    

    BranchList *bl; 
    
    int i,j,k, aktlast;
    
    bl = (BranchList *) activebrunchList;
    
    
    
    
    
    


    aktlast = n;
    

    for( i = 0; i < aktlast ; i++ )
    {
	if (bl[i].minmaxdist < bl[aktlast-1].mindist)
		
		
		
		
	    for( j = 0; (j < aktlast) ; j++ )
		if ((i!=j) && (bl[j].mindist>bl[i].minmaxdist))
		{
		    aktlast = j;

			
			if (j <= i)
				error("Error in pruneBrunchList\n", true);
		    break;
		}
    }

    

    
    
    
    
    


    for( i = 0; i < aktlast; i++)
	if (*nearest_distanz > bl[i].minmaxdist)
	    *nearest_distanz = bl[i].minmaxdist;

    
    
    
    
    
    
    

    for( i = 0; (i < aktlast) && (*nearest_distanz >= bl[i].mindist) ; i++);
    
    aktlast = i;

    

    return (aktlast);
   
}

double objectDIST(double *p1, double *p2)
{

    
    
    
    
    
    
    
    

    double summe = 0;
    int i;
    
    for( i = 0; i < DIMENSION; i++)
	summe += pow(p1[i] - p2[i], 2 );

    
    return(summe);
}

/*****************************************************************
this function returns the maxdist of 2 mbrs
para:
m1: the bounces of the 1st mbr
m2: the bounces of the 2nd mbr
dim: dimensionality
*****************************************************************/

double MbrMAXDIST(double *_m1, double *_m2, int _dim)
{
	double dist=0;
	for (int i=0; i<_dim; i++)
	{
		double d1=fabs(_m1[2*i]-_m2[2*i+1]);
		double d2=fabs(_m1[2*i+1]-_m2[2*i]);
		double d=max(d1, d2);
		dist+=pow(d, 2);
	}
	return dist;
}

/*****************************************************************
this function returns the mindist of 2 mbrs
para:
m1: the bounces of the 1st mbr
m2: the bounces of the 2nd mbr
dim: dimensionality
*****************************************************************/

double MbrMINDIST(double *_m1, double *_m2, int _dim)
{
	double dist=0;
	for (int i=0; i<_dim; i++)
	{
		if (_m1[2*i]>_m2[2*i+1])
			dist+=pow(_m1[2*i]-_m2[2*i+1], 2);
		else if (_m1[2*i+1]<_m2[2*i])
			dist+=pow(_m1[2*i+1]-_m2[2*i], 2);
	}
	return dist;
}

double MINDIST(double *p, double *bounces, int _dim)
{
    
    
    
    

    double summe = 0.0;
    double r;
    int i;

    for(i = 0; i < _dim; i++)
    {
	if (p[i] < bounces[2*i])
	    r = bounces[2*i];
	else
	{
	    if (p[i] > bounces[2*i +1])
		r = bounces[2*i+1];
	    else 
		r = p[i];
	}    

	summe += pow( p[i] - r , 2);
    }
    return(summe);
}

double MAXDIST(double *p, double *bounces, int dim)
{
    double summe = 0.0;
    double r;
    int i;

    for(i = 0; i < dim; i++)
    {
		if (p[i] < bounces[2*i])
			r = bounces[2*i+1];
		else
		{
			if (p[i] > bounces[2*i+1])
			r = bounces[2*i];
			else if (p[i] - bounces[2*i] > bounces[2*i+1] -p[i])
					 r = bounces[2*i];
			     else
					 r = bounces[2*i+1];
		}    

		summe += pow(p[i]-r,2);
    }
    
	return(summe);
}


double MINMAXDIST(double *p, double *bounces)
{

    
    
    
    
    
    
    
    
    
    
    double summe = 0;
    double minimum = 1.0e20;
    double S = 0;

    double rmk, rMi;
    int k,i;

    for( i = 0; i < DIMENSION; i++) 
    { 
	rMi = (	p[i] >= (bounces[2*i]+bounces[2*i+1])/2 )
	    ? bounces[2*i] : bounces[2*i+1];
	S += pow( p[i] - rMi , 2 );
    }

    for( k = 0; k < DIMENSION; k++)
    {  
	
	rmk = ( p[k] <=  (bounces[2*k]+bounces[2*k+1]) / 2 ) ? 
	    bounces[2*k] : bounces[2*k+1];

	summe = pow( p[k] - rmk , 2 );	
	
	rMi = (	p[k] >= (bounces[2*k]+bounces[2*k+1]) / 2 )
	    ? bounces[2*k] : bounces[2*k+1];

	summe += S - pow( p[k] - rMi , 2 );
	
	minimum = min( minimum,summe);
    }

    return(minimum);
    

}
