/*****************************************************************
 * l  i  s  t  .  c  p  p                                        *
 * linked list functions                                         *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/

#include "list.h"


IsectList::IsectList(void) {
  first = NULL;
  current = first;
}

IsectList::~IsectList(void){

	ListItem *temp;

	while (first != NULL && first->next != NULL){
		temp = first->next;		
		delete first;
		first = temp;
	}


	delete first;
	first = NULL;

}


void
IsectList::addItem(Isect *item){
  ListItem *newitem;

  newitem = new ListItem;
  newitem->node = item;
  newitem->next = first;
  first = newitem;

}

bool
IsectList::isEmpty(void){
  if (first==NULL) return true;
  else return false;
}

void
IsectList::beginIteration(void){
  current = first;
}

Isect*
IsectList::getItem(void){
  Isect *temp;

  if(current == NULL) return NULL;
  temp = current->node;
  
  current = current->next;
  return temp;

}
