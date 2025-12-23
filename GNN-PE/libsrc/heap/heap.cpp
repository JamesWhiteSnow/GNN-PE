#include "heap/heap.h"
#include "func/gendef.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>


HeapEntry::HeapEntry()
{
}

HeapEntry::~HeapEntry()
{
}

void HeapEntry::init_HeapEntry(int _dim)
{
	dim = _dim;
}

void HeapEntry::copy(HeapEntry *_he)
{
	key = _he -> key;
	level = _he -> level;
	son1 = _he -> son1;
	son2 = _he -> son2;
}



Heap::Heap()
{
	cont = NULL;
}

Heap::~Heap()
{
	
	delete [] cont;
	cont = NULL;
}

void Heap::enter(HeapEntry *_he, int _pos)

{

	for (int i = used - 1; i >= _pos; i --)
	{
		cont[i + 1].copy(&(cont[i]));
	}
	cont[_pos].copy(_he);
	used ++;

	if (maxused<used)
		maxused=used;
}

void Heap::insert(HeapEntry *_he)
{
	int pos = used;  

	enter(_he, pos);
	
	pos++;
	int parent = pos;
	while (parent != 1)
	{
		int child = parent;
		parent /= 2;
		if (cont[parent - 1].key > cont[child - 1].key)
		{
			HeapEntry *the = new HeapEntry();
			the -> init_HeapEntry(cont[parent - 1].dim);
			the -> copy(&(cont[parent - 1]));
			cont[parent - 1].copy(&(cont[child - 1]));
			cont[child - 1].copy(the);
			delete the;
		}
		else 
			parent = 1;
	}

	if (used > hsize)  
		
		
	{
		error("heap exceeded...\n", true);
	}

	
	
}

bool Heap::remove(HeapEntry *_he)


{
	if (used==0) 
		return false;
	_he -> copy(&(cont[0]));
	used--;
	cont[0].copy(&(cont[used]));
	int parent = 1;
	while (2 * parent <= used)
	{
		int child = 2 * parent;
		if (2 * parent + 1 > used)
			child = 2 * parent;
		else
			if (cont[2 * parent - 1].key < cont[2 * parent].key)
				child = 2 * parent;
			else 
				child = 2 * parent + 1;

		if (cont[parent - 1].key > cont[child - 1].key)
		{
			HeapEntry *the = new HeapEntry();
			the -> init_HeapEntry(cont[parent - 1].dim);
			the -> copy(&(cont[parent - 1]));
			cont[parent - 1].copy(&(cont[child - 1]));
			cont[child - 1].copy(the);
			delete the;
			parent = child; 
		}
		else
			parent = used;
	}

	
	
	return true;
};

void Heap::clean(double _dist)


{
	for (int i = 0; i < used; i ++)
	{
		if (cont[i].key > _dist)
			used = i;
	}
}

/*****************************************************************
this function checks the integrity of the heap. it is an auxiliary
function for debugging

Coded by Yufei Tao 09/01/02
*****************************************************************/

bool Heap::check()
{
	for (int i = 0; i < used; i ++)
	{
		if (cont[i].son1<0 || cont[i].son2<0)
			return false;
	}
	return true;
}

void Heap::init(int _dim, int _hsize)
{
	if (cont)
		delete [] cont;
	hsize = _hsize;
	cont = new HeapEntry [hsize + 1];   
	for (int i = 0; i < hsize + 1; i ++)
		cont[i].init_HeapEntry(_dim);
	used = 0;
	maxused=0;
}

	