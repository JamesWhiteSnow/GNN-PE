#pragma once

#define MAX_HEAP_SIZE 100000

class HeapEntry;
class Heap;

class HeapEntry
{
public:
	int dim;
	int level;
	int son1;
	int son2;
	double key;

	
	HeapEntry();
	~HeapEntry();
	void init_HeapEntry(int _dim);
	void copy(HeapEntry *_he);
};

class Heap
{
public:
	int b;            
	int hsize;        
	int used;         
	int maxused;
	HeapEntry *cont;  

	
	Heap();
	~Heap();
	bool check();
	void clean(double _dist);
	void enter(HeapEntry *_he, int _pos);
	void init(int _dim, int _hsize=MAX_HEAP_SIZE);
	void insert(HeapEntry *_he);
	bool remove(HeapEntry *_he);
};
