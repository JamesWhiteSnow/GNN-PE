#pragma once

#include <stdio.h>
#include "func/gendef.h"


class BlockFile
{
public:
   FILE* fp;			
   char* filename;		 
   int blocklength;	    
   int act_block; 	    
   int number;		    
   bool new_flag;		

   
   BlockFile(char* name, int b_length);
   			       
   ~BlockFile();

   void put_bytes(char* bytes,int num)
		{ fwrite(bytes,num,1,fp); }

   void get_bytes(char* bytes,int num)	     
		{ fread(bytes,num,1,fp); }

   void fwrite_number(int num);	

   int fread_number();		

   void seek_block(int bnum)    
		{ fseek(fp,(bnum-act_block)*blocklength,SEEK_CUR); }

   void read_header(char * header);

   void set_header(char* header);
   					
   bool read_block(Block b,int i);	

   bool write_block(Block b,int i);

   int append_block(Block b);	

   bool delete_last_blocks(int num);

   bool file_new()			
		{ return new_flag; }

   int get_blocklength()	
		{ 
	       return blocklength; }

   int get_num_of_blocks()	
		{ return number; }
};

class CachedBlockFile : public BlockFile
{
public:
   enum uses {free,used,fixed};
   int ptr;		       
   int cachesize;		
   int page_faults;     

   int *cache_cont;	    
   uses *fuf_cont; 		
   int  *LRU_indicator;
   bool  *dirty_indicator;  

   char **cache;   		

   CachedBlockFile(char* name, int blength, int csize);
   					
   ~CachedBlockFile();

   int next();		

   int in_cache(int index);	
   
   bool read_block(Block b,int i);

   bool write_block(Block b,int i);

   bool fix_block(int i);

   bool unfix_block(int i);

   void unfix_all();			

   void set_cachesize(int s);

   void flush();
};
