#pragma once


#include "func/gendef.h"
#include "heap/heap.h"

class LinList;
class SortedLinList;
class Cache;
class RTNode;
class Entry;

class RTree : public Cacheable
{
public:

	int dimension;                       
	int num_of_data;	                 
    int num_of_dnodes;	                 
    int num_of_inodes;	                 
	int root;                            
	bool root_is_data;                   

	RTNode *root_ptr;
    bool *re_level;  
    LinList *re_data_cands; 
	LinList *deletelist;

	int na[10];

	RTree(char *fname, int _b_length, Cache* c, int _dimension);
    RTree(char *fname, Cache* c);
    RTree(char *inpname, char *fname, int _blength, Cache* c, int _dimension);
    ~RTree();

	bool delete_entry(Entry *d);
	bool FindLeaf(Entry *e);
    int get_num() { return num_of_data; }
	void insert(Entry *d);
	void load_root();  
	void NNQuery(double *QueryPoint, SortedLinList *res);
	void rangeQuery(double *mbr, SortedLinList *res);
	void read_header(char *buffer);      
	void write_header(char *buffer);  

	double get_score_linear(Entry *_e, double *_weight, int _dimension);
	double get_score_linear(double *_mbr, double *_weight, int _dimension);
	double get_score_range (Entry *_e, double *_weight, int _dimension, double *_range);
	double get_score_linear(Entry *_e, double *_weight, int _dimension, 
							  double SCORE_FUNC(const double *, const double *, int));
	void rank_qry_constr_linear(double *_weight, double *_qmbr, int _k, Heap *_hp, int *_rslt);
	void rank_qry_inquiry(double *_weight, double _qscore, int *_rslt);
	void rank_qry_linear(double *_weight, int _k, Heap *_hp, int *_rslt);
	void rank_qry_monotone(double *_weight, int _k, Heap *_hp, int *_rslt, 
							double SCORE_FUNC(const double *, const double *, int));

	void kNN(double *_q, int _k, Heap *_hp, double &fardist);
};	
