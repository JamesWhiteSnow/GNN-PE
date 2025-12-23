#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rtree/rtnode.h"
#include "rtree/rtree.h"
#include "rtree/entry.h"
#include "blockfile/blk_file.h"
#include "blockfile/cache.h"
#include "linlist/linlist.h"
#include "func/gendef.h"

RTNode::RTNode(RTree *rt)
  
{
    char *b;
    int header_size;
    Entry * d;
    int i;

    my_tree = rt;
    dimension = rt->dimension;
    num_entries = 0;
	dirty = TRUE;

    d = new Entry();
	d -> init_entry(dimension, NULL);
    header_size = sizeof(char) + sizeof(int);  
    capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();
    delete d;

    entries = new Entry[capacity];
    for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);

	
    b = new char[rt -> file -> get_blocklength()];
    block = rt -> file -> append_block(b);
    delete [] b;
}

RTNode::RTNode(RTree *rt, int _block)
  
{
    char *b;
    int header_size;
    Entry * d;
    int i;

    my_tree = rt;
    dimension = rt->dimension;
    num_entries = 0;
	dirty = FALSE;

    d = new Entry();
	d -> init_entry(dimension, NULL);
    header_size = sizeof(char) + sizeof(int);
    capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();
    delete d;

    entries = new Entry[capacity];
    for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);
    
    block = _block;
    b = new char[rt -> file -> get_blocklength()];
    if (rt -> cache == NULL) 
        rt -> file -> read_block(b, block);
    else
        rt -> cache -> read_block(b, block, rt);

    read_from_buffer(b);
    delete [] b;
}

RTNode::~RTNode()
{
    char *b;

    if (dirty)
    {
        b = new char[my_tree->file->get_blocklength()];
        write_to_buffer(b);

        if (my_tree->cache == NULL) 
            my_tree->file->write_block(b, block);
        else
            my_tree->cache->write_block(b, block, my_tree);

        delete [] b;
    }

    delete [] entries;
}

int RTNode::choose_subtree(double *mbr)
{
    int i, j, follow, minindex, *inside, inside_count, *over;
    double *bmbr, old_o, o, omin, a, amin, f, fmin;

    inside_count = 0;
    inside = new int[num_entries];
    over = new int[num_entries];
    for (i = 0; i < num_entries; i++)
    {
    	switch (entries[i].section(mbr))
    	{
        	case INSIDE:
        	    inside[inside_count++] = i;
        	    break;
        }
    }

    if (inside_count == 1)
        
    	follow = inside[0];
    else if (inside_count > 1)
    
    
    {

		fmin=area(dimension, entries[inside[0]].bounces); minindex=0;
		

    	for (i = 1; i < inside_count; i++)
    	{
    	    f = area(dimension, entries[inside[i]].bounces);
    	    if (f < fmin)
    	    {
    	    	minindex = i;
          		fmin = f;
       	    }
       	}
    	follow = inside[minindex];
    }
    else
    
    
    
    {
       	if (level == 1) 
    	{
            omin = MAXREAL;
    	    fmin = MAXREAL;
    	    amin = MAXREAL;
    	    for (i = 0; i < num_entries; i++)
    	    {
        		enlarge(dimension, &bmbr, mbr, entries[i].bounces);

        		
        		a = area(dimension, entries[i].bounces);
        		f = area(dimension, bmbr) - a;

        		
        		old_o = o = 0.0;
        		for (j = 0; j < num_entries; j++)
        		{
        		    if (j != i)
        		    {
    			        old_o += overlap(dimension,
    					 entries[i].bounces,
    					 entries[j].bounces);
    			        o += overlap(dimension,
    				     bmbr,
    				     entries[j].bounces);
    		        }
    	        }
    	        o -= old_o;

    	        
    	        if ((o < omin) ||
    		    (o == omin && f < fmin) ||
    		    (o == omin && f == fmin && a < amin))
    	        {
    	       	    minindex = i;
        		    omin = o;
        		    fmin = f;
        		    amin = a;
        	    }
    	        delete [] bmbr;
    	    }
        }
        else 
        {
    	    fmin = MAXREAL;
    	    amin = MAXREAL;
    	    for (i = 0; i < num_entries; i++)
    	    {
    	        enlarge(dimension, &bmbr, mbr, entries[i].bounces);

    	        
    	        a = area(dimension, entries[i].bounces);
    	        f = area(dimension, bmbr) - a;

    	        
    	        if ((f < fmin) || (f == fmin && a < amin))
    	        {
    	       	    minindex = i;
    		        fmin = f;
    	            amin = a;
    	        }
	            delete [] bmbr;
	        }
        }

    	follow = minindex;

    	dirty = TRUE;
    }

    delete [] inside;
    delete [] over;

    return follow;
}

R_DELETE RTNode::delete_entry(Entry *e)
{
	RTNode *succ;
	double *tmp;
	if (level > 0)
	{
		if (this == my_tree->root_ptr)
			
		{
			for (int i = 0; i < num_entries; i++)
			{
				tmp = overlapRect(dimension, entries[i].bounces, e -> bounces);
				if (tmp != NULL)
				{
					delete [] tmp;
					succ = entries[i].get_son();
					R_DELETE del_ret;
					del_ret = succ -> delete_entry(e);
					if (del_ret != NOTFOUND)
					{
						switch (del_ret)
						{
						case NORMAL:

							double *mbr;
							mbr = succ -> get_mbr();
							memcpy(entries[i].bounces, mbr, sizeof(double) * 2 * dimension);
							dirty = true;
							delete [] mbr;

							delete entries[i].son_ptr;
							entries[i].son_ptr = NULL;

							return NORMAL;
							break;

						case ERASED:
							delete entries[i].son_ptr;
							entries[i].son_ptr = NULL;

							int j;
							for (j = i; j < num_entries - 1; j++)
								entries[j] = entries[j+1];
							for (j = num_entries - 1; j < capacity; j++)
								entries[j].son_ptr = NULL;

							num_entries--;

							dirty = true;
							return NORMAL;
							break;
						}
					}
				}
			}
			
			return NOTFOUND;
		}
		else
		{
			for (int i = 0; i < num_entries; i++)
			{
				tmp = overlapRect(dimension, entries[i].bounces, e -> bounces);
				if (tmp != NULL)
				{
					delete [] tmp;
					succ = entries[i].get_son();
					R_DELETE del_ret;
					del_ret = succ->delete_entry(e);
					if (del_ret != NOTFOUND)
					{
						switch (del_ret)
						{
						case NORMAL:

							double *mbr;
							mbr = succ -> get_mbr();
							memcpy(entries[i].bounces, mbr, sizeof(double) * 2 * dimension);
							dirty = true;
							delete [] mbr;

							entries[i].del_son();

							return NORMAL;
							break;

						case ERASED:

							entries[i].del_son();

							int j;
							for (j = i; j < num_entries - 1; j++)
								entries[j] = entries[j+1];
							for (j = num_entries - 1; j < capacity; j++)
								entries[j].son_ptr = NULL;
							
							num_entries--;

							dirty = true;
							delete succ;

							if (num_entries < (int)ceil(0.4 * capacity))
							{
								for (int j = 0; j < num_entries; j++)
								{
									Linkable *e;
									e = entries[j].gen_Linkable();
									my_tree -> deletelist -> insert(e);
								}

								my_tree -> num_of_inodes --;
								return ERASED;
							}
							else
								return NORMAL;
							break;
						}
					}
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < num_entries; i++)
		{
			if (entries[i] == (*e))
			{
				my_tree -> num_of_data --;

				for (int j = i; j < num_entries-1; j++)
					entries[j] = entries[j+1];
				
				num_entries--;
				dirty = true;

				if (this != my_tree -> root_ptr && num_entries < (int)ceil(0.4 * capacity))
				{
					for (int k = 0; k < num_entries; k++)
					{
						Linkable *en;
					    en = entries[k].gen_Linkable();
						en -> level = 0;
						my_tree -> deletelist -> insert(en);
					}

					my_tree -> num_of_dnodes --;
					return ERASED;
				}
				else
					return NORMAL;
			}
		}
		return NOTFOUND;
	}
}

void RTNode::enter(Entry *de)
  
{
    if (num_entries > (capacity-1))
        error("RTNode::enter: called, but node is full", TRUE);

    entries[num_entries] = *de;

    num_entries++;

	dirty = true;

    de->son_ptr = NULL;
    delete de;
}

bool RTNode::FindLeaf(Entry *e)
{
	RTNode *succ;
	if (level > 0)
	{
		for (int i = 0; i < num_entries; i++)
		{
			double *f;
			f = overlapRect(my_tree -> dimension,
				  entries[i].bounces, e -> bounces);
			if (f != NULL)
			{
				delete [] f;
				succ = entries[i].get_son();
				bool find;
				find = succ->FindLeaf(e);
				entries[i].del_son();
				if (find)
					return true;
			}
		}
		return false;
	}
	else
	{
		for (int i = 0; i < num_entries; i++)
		{
			if (entries[i] == (*e))
				return true;
		}
		return false;
	}
	return false;
}

double* RTNode::get_mbr()
{
    int i, j;
    double *mbr;

    mbr = new double[2*dimension];
    for (i = 0; i < 2*dimension; i ++ )
        mbr[i] = entries[0].bounces[i];

    for (j = 1; j < num_entries; j++)
    {
    	for (i = 0; i < 2*dimension; i += 2)
    	{
    	    mbr[i]   = min(mbr[i],   entries[j].bounces[i]);
    	    mbr[i+1] = max(mbr[i+1], entries[j].bounces[i+1]);
        }
    }

    return mbr;
}

int RTNode::get_num_of_data()
{
    int i, sum;
    RTNode* succ;

    if (level == 0)
        return num_entries;

    sum = 0;
    for (i = 0; i < num_entries ; i++)
    {
        succ = entries[i].get_son();
        sum += succ->get_num_of_data();
		entries[i].del_son();
    }

    return sum;
}

int RTNode::get_num_of_node()
{
    int i, sum;
    RTNode* succ;

    if (level == 0)
        return 1;

    sum = 1;
    for (i = 0; i < num_entries; i++)
    {
        succ = entries[i].get_son();
        sum += succ->get_num_of_node();
        entries[i].del_son();
    }

    return sum;
}

R_OVERFLOW RTNode::insert(Entry *d, RTNode **sn)
{
    int follow;
    RTNode *succ, *new_succ;
    RTNode *brother;
    Entry *de;
    R_OVERFLOW ret;
    double *mbr,*nmbr;

    int i, last_cand;
    double *center;
    SortMbr *sm;
    Entry *new_entries;

    if (level > 0) 
    {
	  if (level > d -> level)
	  {
        follow = choose_subtree(d -> bounces);

        succ = entries[follow].get_son();

        ret = succ -> insert(d, &new_succ);
    
        mbr = succ -> get_mbr();
        memcpy(entries[follow].bounces, mbr, sizeof(double) * 2 * dimension);
        delete [] mbr;

		entries[follow].del_son();

        if (ret == SPLIT)
        
        {
            if (num_entries == capacity)
         	    error("RTNode::insert: maximum capacity violation", TRUE);

            de = new Entry(dimension, my_tree);
    	    nmbr = new_succ -> get_mbr();
            memcpy(de -> bounces, nmbr, 2 * dimension * sizeof(double));
    	    delete [] nmbr;
            de -> son = new_succ -> block;
			delete new_succ;
            de -> son_ptr = NULL;
            enter(de);

            if (num_entries == (capacity - 1))
            {
        	    brother = new RTNode(my_tree);
        	    my_tree -> num_of_inodes++;
        	    brother -> level = level;
        	    split(brother);
                *sn = brother;
                ret = SPLIT;
        	}
            else
          	    ret = NONE;
        }
        dirty = TRUE;

        return ret;
	  }
	  else 
	  {
		  enter(d);    
		    
		  if (num_entries == (capacity - 1))
            
            
            
		  {
        	brother = new RTNode(my_tree);
        	my_tree -> num_of_inodes++;
        	brother -> level = level;
        	split(brother);
            *sn = brother;
            ret = SPLIT;
		  }
          else
          	ret = NONE;

		  dirty=true;
		  return ret;
	  }	
    }
    else 
    {
        if (num_entries == capacity)
        	error("RTDataNode::insert: maximum capacity violation", TRUE);

        enter(d);

        dirty = TRUE;

        if (num_entries == (capacity - 1))
        
        
        
        {
            if (my_tree->re_level[0] == FALSE && my_tree -> root_ptr -> level != level)
    	    
			
			
			
			
            {
                
                mbr = get_mbr();
                center = new double[dimension];
                for (i = 0; i < dimension; i++)
                     center[i] = (mbr[2*i] + mbr[2*i+1]) / 2.0;

                new_entries = new Entry[capacity];

				for (i = 0; i < capacity; i ++)
					new_entries[i].init_entry(dimension, my_tree);

        	    sm = new SortMbr[num_entries];
        	    for (i = 0; i < num_entries; i++)
        	    {
            		sm[i].index = i;
            		sm[i].dimension = dimension;
            		sm[i].mbr = entries[i].bounces;
            		sm[i].center = center;
                }

                qsort(sm, num_entries, sizeof(SortMbr), sort_center_mbr);

                last_cand = (int) ((double)num_entries * 0.30);

                
                for (i = 0; i < num_entries - last_cand; i++)
    	            new_entries[i] = entries[sm[i].index];

                
                for ( ; i < num_entries; i++)
                {
					Linkable *nd = entries[sm[i].index].gen_Linkable();
                    my_tree -> re_data_cands -> insert(nd);
                }

                
                delete [] entries;
        	    entries = new_entries;
				
        	    delete sm;
        	    delete [] mbr;
        	    delete [] center;
        	    my_tree -> re_level[0] = TRUE;

        	    
        	    num_entries -= last_cand;

        	    
        	    dirty = TRUE;

                return REINSERT;
        	}
           	else  
           	{
        	    *sn = new RTNode(my_tree);
        	    (*sn) -> level = level;
        	    my_tree -> num_of_dnodes++;
        	    split((RTNode *) *sn);
    	    }
    	    return SPLIT;
        }
        else
            return NONE;
    }
}

void RTNode::NNSearch(double *QueryPoint, 
					  SortedLinList *res,
				      double *nearest_distanz)
{
	if (level > 0)
	{
		double *minmax_distanz;		
		int *indexliste;		
		int i,j,k,last,n;
		double akt_min_dist;		
		double minmaxdist,mindist;

		BranchList *activebranchList;
    
		n = num_entries;
    
		
    
		
													
													
		if (res -> get_num() > 0)
		{
			if (*nearest_distanz != res -> get_first() -> distanz)
			{
				printf("testing...\n");
				*nearest_distanz = res -> get_first() -> distanz;
			}
		}

		activebranchList = new BranchList [n]; 
 
		for( i = 0; i < n; i++)
		{
			activebranchList[i].entry_number = i;
			activebranchList[i].minmaxdist = MINMAXDIST(QueryPoint,entries[i].bounces);
			activebranchList[i].mindist = MINDIST(QueryPoint,entries[i].bounces, dimension);	
		}	

		
		qsort(activebranchList,n,sizeof(BranchList),sortmindist); 
 
		
		last = pruneBrunchList(nearest_distanz,activebranchList,n);

		for( i = 0; i < last; i++)
		{
			entries[activebranchList[i].entry_number].get_son()->NNSearch(QueryPoint, res, nearest_distanz);
			entries[i].del_son();
 			
			last = pruneBrunchList(nearest_distanz,activebranchList,last);
		}

		delete [] activebranchList;
	}
	else 
	{
		int i,j;
		double nearest_dist,distanz;
		bool t;
		Linkable *element;

		for (i = 0; i < num_entries; i++)
		{
			
			
			distanz = MINDIST(QueryPoint,entries[i].bounces, dimension);
			
			if (distanz <= *nearest_distanz)
			{
				
				
				
				if (res -> get_num() > 0)
				{
					Linkable *lin = res -> get_first();
					res -> erase();
				}
				
				
				element = entries[i].gen_Linkable();
				element->distanz = distanz;
				res->insert(element);
				
				
				*nearest_distanz = distanz;
			}
		}
	}
}

void RTNode::print()
{
    int i;

	printf("level %d  Block: %d\n", level, block);
	
    for (i = 0; i < num_entries ; i++)
    {
        printf("(%4.1lf, %4.1lf, %4.1lf, %4.1lf)\n",
	       entries[i].bounces[0],
	       entries[i].bounces[1],
	       entries[i].bounces[2],
	       entries[i].bounces[3]);
    }
}

void RTNode::rangeQuery(double *mbr, SortedLinList *res)
{
    int i, n;
    SECTION s;
    RTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
        s = entries[i].section(mbr);
        if (s == INSIDE || s == OVERLAP)
        {
            if (level == 0)
            {
                Linkable *copy;
				copy = entries[i].gen_Linkable();
        		res -> insert(copy);
            }
            else
            {
                succ = entries[i].get_son();
                succ -> rangeQuery(mbr,res);
				entries[i].del_son();
            }
        }
    }
}

void RTNode::read_from_buffer(char *buffer)
{
    int i, j, s;

    
    memcpy(&level, buffer, sizeof(char));
    j = sizeof(char);

    
    memcpy(&num_entries, &(buffer[j]), sizeof(int));
    j += sizeof(int);

    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
    	entries[i].read_from_buffer(&buffer[j]);
    	j += s;
    }
}

int RTNode::split(double **mbr, int **distribution)
{
    bool lu;
    int i, j, k, l, s, n, m1, dist, split_axis;
    SortMbr *sml, *smu;
    double minmarg, marg, minover, mindead, dead, over, *rxmbr, *rymbr;

    n = num_entries;

    m1 = (int) ceil((double)n * 0.40);

    sml = new SortMbr[n];
    smu = new SortMbr[n];
    rxmbr = new double[2*dimension];
    rymbr = new double[2*dimension];

    
    minmarg = MAXREAL;
    for (i = 0; i < dimension; i++)
    
    {
        for (j = 0; j < n; j++)
        {
            sml[j].index = smu[j].index = j;
            sml[j].dimension = smu[j].dimension = i;
            sml[j].mbr = smu[j].mbr = mbr[j];
        }

        
      	qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
        qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

        marg = 0.0;
        
        for (k = 0; k < n - 2 * m1 + 1; k++)
        {
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (l = 0; l < m1 + k; l++)
            {
				for (s = 0; s < 2*dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
				}
			}
			marg += margin(dimension, rxmbr);

			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for ( ; l < n; l++)
            {
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);
        }

        
       	for (k = 0; k < n - 2 * m1 + 1; k++)
        {
            
			
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (l = 0; l < m1+k; l++)
            {
                
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);

            
			
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for ( ; l < n; l++)
            {
                
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);
        }

        if (marg < minmarg)
        {
            split_axis = i;
            minmarg = marg;
        }
    }

    
    for (j = 0; j < n; j++)
    {
		sml[j].index = smu[j].index = j;
		sml[j].dimension = smu[j].dimension = split_axis;
		sml[j].mbr = smu[j].mbr = mbr[j];
    }

    
    qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
    qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

    minover = MAXREAL;
    mindead = MAXREAL;
    
    for (k = 0; k < n - 2 * m1 + 1; k++)
    {
        dead = 0.0;
		for (s = 0; s < 2 * dimension; s += 2)
		{
			rxmbr[s] =    MAXREAL;
			rxmbr[s+1] = -MAXREAL;
		}
		for (l = 0; l < m1 + k; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
				rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
			}
			dead -= area(dimension, sml[l].mbr);
		}
        dead += area(dimension, rxmbr);
		  
		  
		  
		  

		for (s = 0; s < 2*dimension; s += 2)
		{
			rymbr[s] =    MAXREAL;
       		rymbr[s+1] = -MAXREAL;
		}
		for ( ; l < n; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rymbr[s] =   min(rymbr[s],   sml[l].mbr[s]);
				rymbr[s+1] = max(rymbr[s+1], sml[l].mbr[s+1]);
			}
			dead -= area(dimension, sml[l].mbr);
		}
        dead += area(dimension, rymbr);

		over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead)
        {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = TRUE;
        }

		
        dead = 0.0;
		for (s = 0; s < 2*dimension; s += 2)
		{
			rxmbr[s] =    MAXREAL;
			rxmbr[s+1] = -MAXREAL;
		}
		for (l = 0; l < m1+k; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
				rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
			}
			dead -= area(dimension, smu[l].mbr);
		}
        dead += area(dimension, rxmbr);

		for (s = 0; s < 2*dimension; s += 2)
		{
			rymbr[s] =    MAXREAL;
			rymbr[s+1] = -MAXREAL;
		}
		for ( ; l < n; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rymbr[s] =   min(rymbr[s],   smu[l].mbr[s]);
				rymbr[s+1] = max(rymbr[s+1], smu[l].mbr[s+1]);
			}
			dead -= area(dimension, smu[l].mbr);
		}
		
        dead += area(dimension, rymbr);

		over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead)
        {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = FALSE;
        }
    }

    
	
    *distribution = new int[n];
    for (i = 0; i < n; i++)
    {
        if (lu)
            (*distribution)[i] = sml[i].index;
        else
            (*distribution)[i] = smu[i].index;
    }

    delete [] sml;
    delete [] smu;
    delete [] rxmbr;
    delete [] rymbr;

    return dist;
}

void RTNode::split(RTNode *sn)
{
    int i, *distribution, dist, n;
    double **mbr_array;
    Entry *new_entries1, *new_entries2;

    n = num_entries;

    mbr_array = new doubleptr[n];
    for (i = 0; i < n; i++)
       	mbr_array[i] = entries[i].bounces;

    dist = split(mbr_array, &distribution);

    new_entries1 = new Entry[capacity];
    new_entries2 = new Entry[capacity];

	for (i = 0; i < capacity; i ++)
	{
		new_entries1[i].init_entry(dimension, my_tree);
		new_entries2[i].init_entry(dimension, my_tree);
	}

    for (i = 0; i < dist; i++)
       	new_entries1[i] = entries[distribution[i]];

    for (i = dist; i < n; i++)
       	new_entries2[i-dist] = entries[distribution[i]];

    for (i = 0; i < n; i++)
    {
       	entries[i].son_ptr = NULL;
       	sn->entries[i].son_ptr = NULL;
    }
    delete [] entries;
    delete [] sn->entries;

    entries = new_entries1;
    sn->entries = new_entries2;

    num_entries = dist;
    sn->num_entries = n - dist;

    delete [] mbr_array;
	delete [] distribution;
}

void RTNode::write_to_buffer(char *buffer)
{
    int i, j, s;

    
    memcpy(buffer, &level, sizeof(char));
    j = sizeof(char);

    
    memcpy(&buffer[j], &num_entries, sizeof(int));
    j += sizeof(int);

    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
    	entries[i].write_to_buffer(&buffer[j]);
       	j += s;
    }
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

void RTNode::rank_qry_inquiry(double *_weight, double _qscore, int *_rslt)
{
	my_tree->na[level]++;
	if (level==0)
	{
		return; 
	}
	for (int i=0; i<num_entries; i++)
	{
		double range[2];
		my_tree->get_score_range(&(entries[i]), _weight, dimension, range); 
		if (range[0]<=_qscore && _qscore<=range[1])
		{
			RTNode *succ=entries[i].get_son();
			succ->rank_qry_inquiry(_weight, _qscore, _rslt);
			entries[i].del_son();
		}
	}
}

