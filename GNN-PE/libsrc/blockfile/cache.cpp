#include <string.h>
#include "blockfile/cache.h"
#include "blockfile/blk_file.h"

Cache::Cache(int csize, int blength)
{
	int i;

	ptr=0;
	blocklength = blength;

	if (csize >= 0) cachesize=csize;
	else error("Cache size cannot be negative\n", TRUE);

	cache_cont = new int[cachesize];
	cache_tree = new Cacheable*[cachesize];
	fuf_cont = new uses[cachesize];
	LRU_indicator = new int[cachesize];
	dirty_indicator = new bool[cachesize];

	for (i = 0; i < cachesize; i++)
	{
    	cache_cont[i] = 0;
    	cache_tree[i] = NULL;
    	fuf_cont[i] = free;
		LRU_indicator[i] = 0;
		dirty_indicator[i] = false;
	}

	cache = new char*[cachesize];
	for (i = 0; i < cachesize; i++)
		cache[i] = new char[blocklength];

	page_faults = 0;
}

Cache::~Cache()
{
	delete[] cache_cont;
	delete[] fuf_cont;
	delete[] LRU_indicator;
	delete[] cache_tree;
	delete [] dirty_indicator;

	for (int i=0;i<cachesize;i++)
	    delete[] cache[i];
	delete[] cache;
}

int Cache::next()
{
   int ret_val, tmp;

   if (cachesize == 0) return -1;
   else
   {
       if (fuf_cont[ptr] == free)
       {
    	   ret_val = ptr++;  ptr = ptr % cachesize;
    	   return ret_val;
	   }
	   else
	   {
		 tmp= (ptr+1) % cachesize;
		   
		 while (tmp != ptr && fuf_cont[tmp] != free)
		    tmp = (tmp + 1) % cachesize;

		 if (ptr == tmp)	
		 {
	         
             int lru_index = 0; 
             for (int i = 0; i < cachesize; i++)
                if (LRU_indicator[i] > LRU_indicator[lru_index])
                    lru_index = i;

             ptr = lru_index;

    		 if (dirty_indicator[ptr])
				 cache_tree[ptr] -> file -> write_block(cache[ptr], cache_cont[ptr]-1);

    		 fuf_cont[ptr] = free;  dirty_indicator[ptr] = false;

    		 ret_val=ptr++;  ptr = ptr % cachesize;
				
    		 return ret_val;
    	 }
    	 else  
		     return tmp;
		}
   }
}

int Cache::in_cache(int index, Cacheable *rt)
{
   int i;
   int ret_val = -1;
   for (i = 0; i < cachesize; i++)
	   if ((cache_cont[i] == index) && (cache_tree[i] == rt) && (fuf_cont[i] != free))
	   {
	       LRU_indicator[i]=0;
		   
		     
		   ret_val = i;
		     
	   }
	   else if (fuf_cont[i] != free)
           LRU_indicator[i]++;   
   return ret_val;
}

bool Cache::read_block(Block block,int index, Cacheable *rt)
{
	int c_ind;

	index++;	
	if(index <= rt->file->get_num_of_blocks() && index>0)
	{
    	if((c_ind = in_cache(index,rt))>=0)
	   		memcpy(block, cache[c_ind], blocklength);
    	else 
    	{
    	    page_faults++;
    		c_ind = next();
    		if (c_ind >= 0) 
    		{
        		rt -> file -> read_block(cache[c_ind],index-1); 
        		cache_cont[c_ind] = index;
        		cache_tree[c_ind] = rt;
        		fuf_cont[c_ind] = used;
        		LRU_indicator[c_ind] = 0;
        		memcpy(block, cache[c_ind], blocklength);
    		}
    		else
    		    rt -> file -> read_block(block,index - 1); 
    	}
    	return TRUE;
	}
	else
	{
		printf("Requested block %d is illegal.", index - 1);  error("\n", true);
    }

	return false;
}

bool Cache::write_block(Block block, int index, Cacheable *rt)
{
	int c_ind;

	index++;	
	if(index <= rt->file->get_num_of_blocks() && index > 0)
	{
    	c_ind = in_cache(index, rt);
    	if(c_ind >= 0)	
		{
			memcpy(cache[c_ind], block, blocklength);
			dirty_indicator[c_ind] = true;
		}
    	else		
    	{
    		c_ind = next();
    		if (c_ind >= 0)
    		{
        		memcpy(cache[c_ind],block,blocklength);
        		cache_cont[c_ind] = index;
        		cache_tree[c_ind] = rt;
        		fuf_cont[c_ind] = used;
        		LRU_indicator[c_ind] = 0;
				dirty_indicator[c_ind] = true;
    		}
    		else
		  		rt -> file -> write_block(block,index - 1);
	    }
	    return TRUE;
	}
	else
	{
		printf("Requested block %d is illegal.", index - 1);  error("\n", true);
    }

	return false;
}

bool Cache::fix_block(int index, Cacheable *rt)
{
	int c_ind;

	index++;	
	if (index <= rt -> file -> get_num_of_blocks() && index>0)
	{
    	if((c_ind = in_cache(index, rt)) >= 0)
		{
			fuf_cont[c_ind] = fixed;
			return true;
		}
		/*
    	else		
    		if((c_ind=next())>=0)	
    		{
        		rt->file->read_block(cache[c_ind],index-1); 
        		cache_cont[c_ind]=index;
        		cache_tree[c_ind]=rt;
        		fuf_cont[c_ind]=fixed;
    		}
    		else	
    		    return FALSE;
		*/ /*lines commented by TAO Yufei.  this code will read a block from
		    the disk without incrementing the page_faults.  on the second hand,
			we choose not to fetch the page if it is not in memory*/
		else
			return false;
	}
	else
	{
		printf("Requested block %d is illegal.", index - 1);  error("\n", true);
    }

	return false;;
}

bool Cache::unfix_block(int index, Cacheable *rt)
{
	int i;

	i = 0;
	index++;	
	if(index <= rt -> file -> get_num_of_blocks() && index > 0)
	{
    	while(i < cachesize && (cache_cont[i] != index || fuf_cont[i] == free))
    		i++;
    	if (i != cachesize)
    		fuf_cont[i] = used;
    	return TRUE;
	}
	else
	{
		printf("Requested block %d is illegal.", index - 1);  error("\n", true);
    }

	return false;
}

void Cache::unfix_all()
{
	int i;

	for (i = 0; i < cachesize; i++)
		if (fuf_cont[i] == fixed)
			fuf_cont[i] = used;
}

void Cache::set_cachesize(int size)
{
	int i;

	if (size >= 0)
	{
    	ptr = 0;
    	flush();
    	for(i = 0; i < cachesize; i++)
    		delete[] cache[i];
    	delete[] cache;

    	delete[] cache_cont;
    	delete[] cache_tree;
    	delete[] fuf_cont;
    	delete[] LRU_indicator;
		delete[] dirty_indicator;

    	cachesize = size;
    	cache_cont = new int[cachesize];
    	cache_tree = new Cacheable*[cachesize];
    	LRU_indicator = new int[cachesize];
    	fuf_cont = new uses[cachesize];
		dirty_indicator = new bool[cachesize];

    	for (i=0; i<cachesize; i++)
    	{
    		cache_cont[i] = 0;
    		cache_tree[i] = NULL;
    		fuf_cont[i] = free;
			LRU_indicator[i] = 0;
			dirty_indicator[i] = 0;
    	}

    	cache = new char*[cachesize];
    	for (i = 0; i < cachesize; i++)
    		cache[i] = new char[blocklength];
	}
	else
	{
	    error("Cache size cannot be negative\n",TRUE);
	}
}

void Cache::flush()
{
	int i;

	for (i=0; i<cachesize; i++)
	{
		if (fuf_cont[i] != free && dirty_indicator[i])
				cache_tree[i]->file->write_block(cache[i], cache_cont[i]-1); 
		fuf_cont[i] = free;
	}
}