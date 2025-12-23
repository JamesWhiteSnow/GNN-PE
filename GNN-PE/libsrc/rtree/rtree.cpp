#include <math.h>
#include <string.h>
#include "rtree/rtree.h"
#include "rtree/entry.h"
#include "rtree/rtnode.h"
#include "blockfile/cache.h"
#include "blockfile/blk_file.h"
#include "linlist/linlist.h"
#include "heap/heap.h"

RTree::RTree(char *fname, int _b_length, Cache *c, int _dimension)
  
{
    file = new BlockFile(fname, _b_length);
    cache = c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
	  
	  
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr -> block;
}

RTree::RTree(char *fname, Cache *c)
  
{
    int j;

    file = new BlockFile(fname, 0);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    file -> read_header(header);
    read_header(header);
	delete [] header;

    root_ptr = NULL;
}

RTree::RTree(char *inpname, char *fname, int _b_length, Cache *c, int _dimension)
  
{
    Entry *d;
    FILE *fp;
    file = new BlockFile(fname, _b_length);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    read_header(header);
	delete [] header;

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTNode(this);
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr->block;

	int record_count = 0;

    if((fp = fopen(inpname,"r")) == NULL)
    {
      delete this;
      error("Cannot open R-Tree text file", TRUE);
    }
    else
    {
	  int id=0;
	  double x0,y0,x1,y1;
      while (!feof(fp))
      {
		record_count ++;
		





		if (record_count%100 == 0)
		{
			for (int i = 0; i < 79; i ++)  
			  printf("\b");
			printf("inserting object %d", record_count);
		}

		fscanf(fp, "%d ", &id);
		d = new Entry(dimension, NULL);
    	d -> son = id;
		for (int i=0; i<dimension; i++)
		{
    		fscanf(fp, "%f %f ", &(d->bounces[2*i]), &(d->bounces[2*i+1]));
		}

    	insert(d);
		  
      }
    }

	fclose(fp);

	printf("\n");
	delete root_ptr;
	root_ptr = NULL;
}

RTree::~RTree()
{
    int j;

	char *header = new char[file -> get_blocklength()];
    write_header(header);
    file->set_header(header);
    delete [] header;

    if (root_ptr != NULL)
    {
        delete root_ptr;
        root_ptr = NULL;
    }

	if (cache)
      cache -> flush();

    delete file;

    delete re_data_cands;
	delete deletelist;

    printf("This R-Tree contains %d internal, %d data nodes and %d data\n",
	   num_of_inodes, num_of_dnodes, num_of_data);
}

bool RTree::delete_entry(Entry *d)
{
	load_root();

	R_DELETE del_ret;
	del_ret=root_ptr->delete_entry(d);

	if (del_ret == NOTFOUND) return false;
	if (del_ret == ERASED) 
		error("RTree::delete_entry--The root has been deleted\n",true);
 
	if (root_ptr -> level > 0 && root_ptr -> num_entries == 1)
		
		
	{
		root = root_ptr -> entries[0].son;
		delete root_ptr;
		root_ptr = NULL;
		load_root();
		num_of_inodes--;
	}

	
	while (deletelist -> get_num() > 0)
	{
		Linkable *e;
		e = deletelist -> get_first();
		Entry *new_e = new Entry(dimension, NULL);
		new_e -> set_from_Linkable(e);
		deletelist -> erase();
		insert(new_e);
	}

	delete root_ptr;
	root_ptr = NULL;

	return true;
}

bool RTree::FindLeaf(Entry *e)
{
	load_root();
	return root_ptr -> FindLeaf(e);
}

void RTree::insert(Entry* d)
{
    int i, j;
    RTNode *sn;
    RTNode *nroot_ptr;
    int nroot;
    Entry *de;
    R_OVERFLOW split_root;
    Entry *dc;
    double *nmbr;

    
    load_root();

    
    re_level = new bool[root_ptr -> level + 1];
    for (i = 0; i <= root_ptr -> level; i++)
        re_level[i] = FALSE;

    
    
    Linkable *new_link;
	new_link = d -> gen_Linkable();
	re_data_cands -> insert(new_link);

	delete d;  

    j = -1;
    while (re_data_cands -> get_num() > 0)
    {
        
	    Linkable *d_cand;
		d_cand = re_data_cands -> get_first();
        if (d_cand != NULL)
        {
            
            
            
			dc = new Entry(dimension, NULL);
            dc -> set_from_Linkable(d_cand);
            re_data_cands -> erase();

            
			split_root = root_ptr -> insert(dc, &sn);
        }
        else
	        error("RTree::insert: inconsistent list re_data_cands", TRUE);

    	if (split_root == SPLIT)
    	
    	{
    	    nroot_ptr = new RTNode(this);
    	    nroot_ptr -> level = root_ptr -> level + 1;
    	    num_of_inodes++;
    	    nroot = nroot_ptr -> block;

    	    de = new Entry(dimension, this);
    	    nmbr = root_ptr -> get_mbr();
    	    memcpy(de->bounces, nmbr, 2*dimension*sizeof(double));
    	    delete [] nmbr;
    	    de->son = root_ptr->block;
    	    de->son_ptr = root_ptr;
    	    nroot_ptr -> enter(de);

    	    de = new Entry(dimension, this);
    	    nmbr = sn -> get_mbr();
    	    memcpy(de -> bounces, nmbr, 2*dimension*sizeof(double));
    	    delete [] nmbr;
    	    de -> son = sn -> block;
    	    de -> son_ptr = sn;
    	    nroot_ptr->enter(de);

    	    root = nroot;
            root_ptr = nroot_ptr;

            root_is_data = FALSE;
        }
        j++;
    }

    num_of_data++;

    delete [] re_level;

	delete root_ptr;
	root_ptr = NULL;
}

void RTree::load_root()
{
    if (root_ptr == NULL)
        root_ptr = new RTNode(this, root);
}

void RTree::NNQuery(double *QueryPoint,
					SortedLinList *res)
{
      double nearest_distanz;

      
      load_root();

      nearest_distanz = MAXREAL;

      root_ptr->NNSearch(QueryPoint,res,&nearest_distanz);

	  delete root_ptr;
	  root_ptr = NULL;
}

void RTree::rangeQuery(double *mbr, SortedLinList *res)
{
    load_root();

    root_ptr -> rangeQuery(mbr,res);

	delete root_ptr;
	root_ptr = NULL;
}

void RTree::read_header(char *buffer)
{
    int i;

    memcpy(&dimension, buffer, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&num_of_data, &buffer[i], sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&num_of_dnodes, &buffer[i], sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root_is_data, &buffer[i], sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);
}

void RTree::write_header(char *buffer)
{
    int i;

    memcpy(buffer, &dimension, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);
}

/*****************************************************************
this function gets the score of an entry with respect to a linear
preference function
para:
e: the leaf entry
weight: the prefrence vector
dimension: dimensionality
*****************************************************************/

double RTree::get_score_linear(Entry *_e, double *_weight, int _dimension)
{
	double max_score=-MAXREAL; 

	  

	
	double end_cnt=1;
	for (int i=0; i<_dimension; i++)
		end_cnt*=2;
	
	for (int i=0; i<end_cnt; i++)
	{
		
		int modseed=1;
		double score=0;
		for (int j=0; j<_dimension; j++)
		{
			int sub=i/modseed%2;

			score+=_e->bounces[2*j+sub]*_weight[j];

			modseed*=2;
		}
		
		if (score>max_score)
			max_score=score;
		
	}
	


	return max_score;
}

double RTree::get_score_linear(double *_mbr, double *_weight, int _dimension)
{
	double max_score=-MAXREAL; 

	  

	
	double end_cnt=1;
	for (int i=0; i<_dimension; i++)
		end_cnt*=2;
	
	for (int i=0; i<end_cnt; i++)
	{
		
		int modseed=1;
		double score=0;
		for (int j=0; j<_dimension; j++)
		{
			int sub=i/modseed%2;

			score+=_mbr[2*j+sub]*_weight[j];

			modseed*=2;
		}
		
		if (score>max_score)
			max_score=score;
		
	}
	


	return max_score;
}

double RTree::get_score_linear(Entry *_e, double *_weight, int _dimension, 
							  double SCORE_FUNC(const double *, const double *, int))
{
	double max_score=-MAXREAL; 
	double *tuple=new double[dimension]; 
	  

	
	double end_cnt=1;
	for (int i=0; i<_dimension; i++)
		end_cnt*=2;
	
	for (int i=0; i<end_cnt; i++)
	{
		
		int modseed=1;
		double score=0;
		for (int j=0; j<_dimension; j++)
		{
			int sub=i/modseed%2;
			tuple[j]=_e->bounces[2*j+sub];

			modseed*=2;
		}
		score=SCORE_FUNC(tuple, _weight, _dimension);
		
		if (score>max_score)
			max_score=score;
		
	}
	

	delete [] tuple;
	return max_score;
}

/*****************************************************************
this function gets the score range of an entry with respect to a linear
preference function
para:
e: the leaf entry
weight: the prefrence vector
dimension: dimensionality
range: the range returned
*****************************************************************/

double RTree::get_score_range (Entry *_e, double *_weight, int _dimension, double *_range)
{
	double max_score=-MAXREAL; 
	double min_score=MAXREAL;

	  

	
	double end_cnt=1;
	for (int i=0; i<_dimension; i++)
		end_cnt*=2;
	
	for (int i=0; i<end_cnt; i++)
	{
		
		int modseed=1;
		double score=0;
		for (int j=0; j<_dimension; j++)
		{
			int sub=i/modseed%2;

			score+=_e->bounces[2*j+sub]*_weight[j];

			modseed*=2;
		}
		
		if (score>max_score)
			max_score=score;
		if (score<min_score)
			min_score=score;
		
	}
	


	_range[0]=min_score; _range[1]=max_score;
	return max_score;
}

/*****************************************************************
this function performs a rank-inquiry query with linear function using
the breadth-first search.
para:
weight: an array (with size equal to the dimensionality) storing
  the query vector
qscore: the query score
rslt: the link list storing the ids of the returned tuples (disabled
for the time being)

Coded by Yufei Tao 05/04/02
*****************************************************************/

void RTree::rank_qry_inquiry(double *_weight, double _qscore, int *_rslt)
{

for (int i=0; i<10; i++)
na[i]=0;

	load_root();
	root_ptr->rank_qry_inquiry(_weight, _qscore, _rslt);
	delete root_ptr;
	root_ptr=NULL;
}

/*****************************************************************
this function performs a constrain ranked query with linear function using
the breadth-first search.
para:
weight: an array (with size equal to the dimensionality) storing
  the query vector
qmbr: the query mbr
k: top "k" query
hp: the heap for the breadth-first search
rslt: the link list storing the ids of the returned tuples

Coded by Yufei Tao 05/04/02
*****************************************************************/

void RTree::rank_qry_constr_linear(double *_weight, double *_qmbr, int _k, Heap *_hp, int *_rslt)
{
	int found_cnt=0;


for (int i=0; i<10; i++)
na[i]=0;

	
	HeapEntry *he=new HeapEntry();
	he->son1=root;
	he->key=0;
	he->level=1; 
	_hp->insert(he);
	delete he;
	

	while(_hp->used>0)
	{
		
		HeapEntry *he=new HeapEntry();
		_hp->remove(he);
		int son=he->son1;
		int level=he->level;
		delete he;
		
		
		if (level==0)
		{
			_rslt[found_cnt]=son;
			found_cnt++;
			if (found_cnt==_k)
				return;
		}
		else  
		{
			RTNode *child=new RTNode(this, son);
			for (int i=0; i<child->num_entries; i++)
			{
				double *ovrp=overlapRect(dimension, child->entries[i].bounces, _qmbr);
				if (ovrp)
				{
					
					
					double score=get_score_linear(ovrp, _weight, dimension);
					delete []ovrp;
					
					HeapEntry *he=new HeapEntry();
					he->son1=child->entries[i].son;
					he->level=child->level;
					he->key=-score; 
					  
					  
					  
					_hp->insert(he);
					delete he;
					
				}
			}
na[child->level]++;
			delete child;
		}
	}
}

/*****************************************************************
this function performs a ranked query with linear function using
the breadth-first search.
para:
weight: an array (with size equal to the dimensionality) storing
  the query vector
k: top "k" query
hp: the heap for the breadth-first search
rslt: the link list storing the ids of the returned tuples

Coded by Yufei Tao 20/03/02
*****************************************************************/

void RTree::rank_qry_linear(double *_weight, int _k, Heap *_hp, int *_rslt)
{
	int found_cnt=0;


for (int i=0; i<10; i++)
na[i]=0;

	
	HeapEntry *he=new HeapEntry();
	he->son1=root;
	he->key=0;
	he->level=1; 
	_hp->insert(he);
	delete he;
	

	while(_hp->used>0)
	{
		
		HeapEntry *he=new HeapEntry();
		_hp->remove(he);
		int son=he->son1;
		int level=he->level;
		delete he;
		
		
		if (level==0)
		{
			_rslt[found_cnt]=son;
			found_cnt++;
			if (found_cnt==_k)
				return;
		}
		else
		{
			RTNode *child=new RTNode(this, son);
			for (int i=0; i<child->num_entries; i++)
			{
				
				double score=get_score_linear(&(child->entries[i]), _weight, dimension);
				
				HeapEntry *he=new HeapEntry();
				he->son1=child->entries[i].son;
				he->level=child->level;
				he->key=-score; 
				  
				  
				  
				_hp->insert(he);
				delete he;
				
			}
na[child->level]++;
			delete child;
		}
	}
}

/*****************************************************************
this function performs a ranked query with a monotone function using
the breadth-first search.
para:
weight: an array (with size equal to the dimensionality) storing
  the query vector
k: top "k" query
hp: the heap for the breadth-first search
rslt: the link list storing the ids of the returned tuples
SCORE_FUNC: the score function used to evaluate a tuple

Coded by Yufei Tao 20/03/02
*****************************************************************/

void RTree::rank_qry_monotone(double *_weight, int _k, Heap *_hp, int *_rslt, 
							double SCORE_FUNC(const double *, const double *, int))
{
	int found_cnt=0;


for (int i=0; i<10; i++)
na[i]=0;

	
	HeapEntry *he=new HeapEntry();
	he->son1=root;
	he->key=0;
	he->level=1; 
	_hp->insert(he);
	delete he;
	

	while(_hp->used>0)
	{
		
		HeapEntry *he=new HeapEntry();
		_hp->remove(he);
		int son=he->son1;
		int level=he->level;
		delete he;
		
		
		if (level==0)
		{
			_rslt[found_cnt]=son;
			found_cnt++;
			if (found_cnt==_k)
				return;
		}
		else
		{
			RTNode *child=new RTNode(this, son);
			for (int i=0; i<child->num_entries; i++)
			{
				
				double score=get_score_linear(&(child->entries[i]), _weight, dimension, SCORE_FUNC);
				
				HeapEntry *he=new HeapEntry();
				he->son1=child->entries[i].son;
				he->level=child->level;
				he->key=-score; 
				  
				  
				  
				_hp->insert(he);
				delete he;
				
			}
na[child->level]++;
			delete child;
		}
	}
}

/*****************************************************************
this function performs a kNN search using the best-first algo.
para:
q: the query point
k: top "k" query
hp: the heap for the best-first search
fardist: the farthest distance from the query to the k-th NN

Coded by Yufei Tao 16/12/02
*****************************************************************/

void RTree::kNN(double *_q, int _k, Heap *_hp, double &fardist)
{
	fardist=MAXREAL;
	int found_cnt=0;


for (int i=0; i<10; i++)
na[i]=0;

	
	
	HeapEntry *he=new HeapEntry();
	he->son1=root;
	he->key=0;
	he->level=1; 
	_hp->insert(he);
	delete he;
	

	while(_hp->used>0)
	{
		
		HeapEntry *he=new HeapEntry();
		_hp->remove(he);
		int son=he->son1;
		int level=he->level;
		if (level==0)
		fardist=he->key;
		delete he;
		
		
		if (level==0)
		{

			found_cnt++;
			if (found_cnt==_k)
				return;
		}
		else
		{
			RTNode *child=new RTNode(this, son);
			for (int i=0; i<child->num_entries; i++)
			{
				
				HeapEntry *he=new HeapEntry();
				he->son1=child->entries[i].son;
				he->level=child->level-1;
				he->key=MINDIST(_q, child->entries[i].bounces, dimension); 
				_hp->insert(he);
				delete he;
				
			}
na[child->level]++;
			delete child;
		}
	}
}



