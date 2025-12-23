#pragma once


#include "func/gendef.h"

class Cacheable;

class Cache
{
public:
   enum uses {free,used,fixed};	
   int ptr;		        
   int cachesize;	
   int blocklength;
   int page_faults;
   int *cache_cont;	 
   Cacheable **cache_tree;  
   uses *fuf_cont; 		
   int  *LRU_indicator; 
   bool *dirty_indicator;  
	   
   char **cache;   		

   
   int next();		

   int in_cache(int index, Cacheable *rt);

   Cache(int csize, int blength);

   ~Cache();

   bool read_block(Block b, int i, Cacheable *rt);

   bool write_block(Block b, int i, Cacheable *rt);

   bool fix_block(int i, Cacheable *rt);

   bool unfix_block(int i, Cacheable *rt);

   void unfix_all();

   void set_cachesize(int s);

   void flush();			
};

