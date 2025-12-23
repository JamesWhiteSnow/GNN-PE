#pragma once


#include "func/gendef.h"

class RTNode;
class RTree;
struct Linkable;

class Entry 
{
public:

	int son;
	double *bounces;                     

	int dimension;                      
	int level;
    RTree *my_tree;                     
    RTNode *son_ptr;                    
   


	Entry();
	Entry(int dimension, RTree *rt);
    ~Entry();

	void del_son();
	Linkable *gen_Linkable();
	int get_size(); 
	RTNode *get_son();
	void init_entry(int _dimension, RTree *_rt);
	void read_from_buffer(char *buffer);
    SECTION section(double *mbr);        
	bool section_circle(double *center, double radius);
	void set_from_Linkable(Linkable *link);
    void write_to_buffer(char *buffer); 

    virtual Entry & operator = (Entry &_d);
	bool operator == (Entry &_d);
};