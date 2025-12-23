#pragma once


#include "func/gendef.h"

class SortedLinList;
class Entry;
class RTree;
class Heap;

class RTNode
{
public:

	char level; 
	int block;
	int num_entries;
	Entry *entries;

	bool dirty;
	int capacity;
    int dimension;
	RTree *my_tree;  

	RTNode(RTree *rt);
    RTNode(RTree *rt, int _block);
    ~RTNode();

    int choose_subtree(double *brm);
	R_DELETE delete_entry(Entry *e); 
	void enter(Entry *de);
	bool FindLeaf(Entry *e);
	double *get_mbr();
	int get_num_of_data();
	int get_num_of_node();
	R_OVERFLOW insert(Entry *d, RTNode **sn);
	bool is_data_node() { return (level==0) ? TRUE : FALSE ;};
	void NNSearch(double *QueryPoint, SortedLinList *res,
				      double *nearest_distanz);
	void print();
	void rangeQuery(double *mbr, SortedLinList *res);
    void read_from_buffer(char *buffer);
	int split(double **mbr, int **distribution);
	void split(RTNode *sn);
	void write_to_buffer(char *buffer); 
	
	void rank_qry_inquiry(double *_weight, double _qscore, int *_rslt);
};