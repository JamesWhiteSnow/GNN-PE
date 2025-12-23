#include "linlist/linlist.h"
#include "func/gendef.h"

SLink::SLink()
{
    d = NULL;
    next = prev = NULL;
}


SLink::~SLink()
{
    delete d;
}


LinList::LinList()
{
    anz = 0;
    akt_index = -1;
    akt = first = last = NULL;
}


LinList::~LinList()
{
    SLink *hf;

    
    akt = first;
    while (akt != NULL)
    {
	hf = akt->next;
	delete akt;
	akt = hf;
    }
}


void LinList::insert(Linkable *f)
{
    SLink *sd;

    
    sd = new SLink;
    sd->d = f;

    
    sd->next = first;
    sd->prev = NULL;
    if (first != NULL)
	first->prev = sd;

    
    anz++;
    first = sd;
    if (last == NULL)
	last = sd;

    
    akt = NULL;
    akt_index = -1;
}


bool LinList::erase()
{
    SLink *n_akt;

    
    if (akt)
    {
	
	if (akt == first)
	{
	    
	    if (akt == last)
	    {
		akt_index = -1;
		first = last = NULL;
		n_akt = NULL;
	    }
	    else
	    {
		
		(akt->next)->prev = NULL;
		first = akt->next;
		n_akt = first;
		akt_index = 0;
	    }
	}
	else
	{
	    
	    if (akt == last)
	    {
		(akt->prev)->next = NULL;
		last = akt->prev;
		n_akt = NULL;
		akt_index = -1;
	    }
	    else
	    
	    {
		(akt->next)->prev = akt->prev;
		(akt->prev)->next = akt->next;
		n_akt = akt->next;
		akt_index++;
	    }
	}

	
	delete akt;

	
	akt = n_akt;

	
	anz--;
	return TRUE;
    }
    return FALSE;
}


Linkable* LinList::get(int i)

{
    bool ahead;   
    int j;

    
    if (i >= anz)
	return NULL;

    
    if (i == akt_index)
	return akt->d;

    
    if (akt_index == -1)
    {
	
	if (i < (anz / 2))
	{
	    akt = first;
	    akt_index = 0;
	    ahead = TRUE;
	}
	else
	{
	    akt = last;
	    akt_index = anz - 1;
	    ahead = FALSE;
	}
    }
    else
    {
	
	if (i < akt_index)
	{
	    
	    if ((akt_index - i) > i)
	    {
		akt = first;
		akt_index = 0;
		ahead = TRUE;
	    }
	    else
		ahead = FALSE;
	}
	else
	{
	    
	    if ((i - akt_index) > ((anz-1) - i))
	    {
		akt = last;
		akt_index = anz - 1;
		ahead = FALSE;
	    }
	    else
		ahead = TRUE;
	}
    }


    
    if (ahead)
    {
	for (j = akt_index; j < i; j++)
	{
	    if (!akt)
		error("LinList::get: List seems to be inkonsistent", TRUE);
	    akt = akt->next;
	}
    }
    else
    {
	for (j = akt_index; j > i; j--)
	{
	    if (!akt)
		error("LinList::get: List seems to be inkonsistent", TRUE);
	    akt = akt->prev;
	}
    }

    akt_index = i;
    return akt->d;
}


Linkable* LinList::get_first()
{
    akt = first;

    if (akt != NULL)
    {
	akt_index = 0;
	return akt->d;
    }
    else
	return NULL;
}


Linkable* LinList::get_last()
{
    akt = last;

    if (akt != NULL)
    {
	akt_index = anz - 1;
	return akt->d;
    }
    else
	return NULL;
}


Linkable* LinList::get_next()
{
    akt = akt->next;

    if (akt != NULL)
    {
	akt_index++;
	return akt->d;
    }
    else
    {
	akt_index = -1;
	return NULL;
    }
}


Linkable* LinList::get_prev()
{
    akt = akt->prev;

    if (akt != NULL)
    {
	akt_index--;
	return akt->d;
    }
    else
    {
	akt_index = -1;
	return NULL;
    }
}

void LinList::print()
{
    SLink *cur = first;

    while (cur != NULL)
    {
	  printf("%d %f %f %f %f\n", cur->d->son, cur->d->bounces[0], cur->d->bounces[1], cur->d->bounces[2],  cur->d->bounces[3]);
	  cur = cur->next;
    }
}

void LinList::check()
{
    SLink *f, *old_f;
    int myanz;
    char buffer[255];

    old_f = first;
    
    if (old_f == NULL)
    {
	if (last != NULL)
	    error("LinList::check: first == NULL, last != NULL", FALSE);
	if (anz != 0)
	    error("LinList::check: first == NULL, anz != 0", FALSE);
	return;
    }

    myanz = 1;
    if (old_f->prev != NULL)
    {
	error("LinList::check: Listenkopf.prev ungleich NULL", FALSE);
	return;
    }

    for (f = old_f->next; f != NULL; f = f->next)
    {
	if (f->prev != old_f)
	{
	    error("LinList::check: Rueckwaertsverkettung fehlerhaft", FALSE);
            return;
        }
	if (old_f->next != f)
	{
	    error("LinList::check: Vorwaertsverkettung fehlerhaft", FALSE);
            return;
        }
	old_f = f;

	myanz ++;
	if (myanz > anz)
	{
	    sprintf(buffer, "LinList::check: anz (%d != %d) passt nicht", myanz, anz);
	    error(buffer, FALSE);
	    return;
	}
    }

    if (old_f->next != NULL)
    {
	error("LinList::check: Listenende.next ungleich NULL", FALSE);
	return;
    }

    if (last != old_f)
    {
	error("LinList::check: last ungleich Listenende", FALSE);
	return;
    }

    if (myanz != anz)
    {
	sprintf(buffer, "LinList::check: anz (%d != %d) passt nicht", myanz, anz);
	error(buffer, FALSE);
    }
}






SortedLinList::SortedLinList()
    : LinList()
{
    increasing = TRUE;
}



void SortedLinList::set_sorting(bool _increasing)
{
    increasing = _increasing;
}

void SortedLinList::insert(Linkable *f)
{
    SLink *sd;

    
    sd = new SLink;
    sd->d = f;

    
    if (first == NULL)
    {
	first = sd;
	last = sd;
	anz = 1;
	sd->next = sd->prev = NULL;

	return;
    }

    
    if (increasing)
	for (akt = first; akt != NULL && (akt->d->son) < f->son;
				 akt = akt->next)
	    ;
    else
	for (akt = first; akt != NULL && akt->d->son > f->son; akt = akt->next)
	    ;

    
    if (akt != NULL)
    {
	if (akt == first)
	    first = sd;
	else
	    (akt->prev)->next = sd;
	sd->next = akt;
	sd->prev = akt->prev;
	akt->prev = sd;
    }
    else
    
    
    {
	sd->prev = last;
	sd->next = NULL;
	last->next = sd;
	last = sd;
    }

    
    anz ++;
    akt_index = -1;
}



void SortedLinList::sort(bool _increasing)

{
    SLink *old_akt;
    Linkable *s;
    bool any, swap;

    set_sorting(_increasing);

    any = TRUE;
    while (any)
    {
	any = FALSE;
	old_akt = NULL;
	for (akt = first; akt != NULL; akt = akt->next)
	{
	    
	    if (old_akt != NULL)
	    {
		swap = FALSE;
		if (_increasing)
		{
		    if ((akt->d->son) > (old_akt->d->son))
			swap = TRUE;
		}
		else
		{
		    if ((akt->d->son) < (old_akt->d->son))
			swap = TRUE;
		}

		if (swap)
		{
		    
		    s = akt->d;
		    akt->d = old_akt->d;
		    old_akt->d = s;

		    
		    
		    any = TRUE;
		}
	    }
	    old_akt = akt;
	}
    }

    akt_index = -1;
}

