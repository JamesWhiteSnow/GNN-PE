#pragma once


#include <stdio.h>
#include <ctype.h>

#define BFHEAD_LENGTH (sizeof(int)*2)    

#define TRUE 1
#define FALSE 0

#define SEEK_CUR 1
#define SEEK_SET 0
#define SEEK_END 2

typedef char Block[];

#define MAXREAL         1e20
#define doubleZERO       1e-20
#define MAX_DIMENSION   256

#define DIMENSION 2

#define TRUE 1
#define FALSE 0

#define min(a, b) (((a) < (b))? (a) : (b)  )
#define max(a, b) (((a) > (b))? (a) : (b)  )

class BlockFile;  
class Cache;
class Cacheable   
                  
{
public:
	BlockFile *file;
	Cache *cache;
};
  
class CmdIntrpr  
                  
{
public:
	int cnfrm_cmd(char *_msg)
	{ char c = ' ';
	  while (c != 'N' && c != 'Y')
	  { printf("%s (y/n)?", _msg);
	    c = getchar(); 
		char tmp;
		while ((tmp = getchar()) != '\n');
		c = toupper(c); }
	  if (c == 'N') return 0; else return 1; }
  
	void get_cmd(char *_msg, char *_cmd)
	{ printf("%s", _msg);  
	  char *c = _cmd;
	  while (((*c) = getchar()) != '\n')
	    c++;
	  *c = '\0'; } 

	virtual bool build_tree(char *_tree_fname, char *_data_fname, int _b_len, int _dim, int _csize) = 0;
	virtual void free_tree() = 0;
	virtual int qry_sngle(double *_mbr, int *_io_count) = 0;
	/*
	virtual int qry_wrkld(char *_wrkld_fname, int *_io_count, bool _display) = 0;
	*/
	virtual void run() = 0;
	virtual void version() = 0;
};
  
enum SECTION {OVERLAP, INSIDE, S_NONE};
enum R_OVERFLOW {SPLIT, REINSERT, NONE};
enum R_DELETE {NOTFOUND,NORMAL,ERASED};
typedef double *doubleptr;
  
struct SortMbr
{
    int dimension;
    double *mbr;
    double *center;
    int index;
};

struct BranchList
{
    int entry_number;
    double mindist;
    double minmaxdist;
    bool section;
};


void error(char *_errmsg, bool _terminate);

double area(int dimension, double *mbr);
double margin(int dimension, double *mbr);
double overlap(int dimension, double *r1, double *r2);
double* overlapRect(int dimension, double *r1, double *r2);
double objectDIST(double *p1, double *p2);
double MINMAXDIST(double *Point, double *bounces);
double MINDIST(double *Point, double *bounces, int Pdim);
double MAXDIST(double *p, double *bounces, int dim);
double MbrMINDIST(double *_m1, double *_m2, int _dim);
double MbrMAXDIST(double *_m1, double *_m2, int _dim);

bool inside(double &p, double &lb, double &ub);
void enlarge(int dimension, double **mbr, double *r1, double *r2);
bool is_inside(int dimension, double *p, double *mbr);
int pruneBrunchList(double *nearest_distanz, const void *activebrunchList, 
		    int n);
bool section(int dimension, double *mbr1, double *mbr2);
bool section_c(int dimension, double *mbr1, double *center, double radius);

int sort_lower_mbr(const void *d1, const void *d2);
int sort_upper_mbr(const void *d1, const void *d2);
int sort_center_mbr(const void *d1, const void *d2);
int sortmindist(const void *element1, const void *element2);

#ifdef UNIX
void strupr(char *_msg);
#endif
