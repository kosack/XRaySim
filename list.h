/*****************************************************************
 * l  i  s  t  .  h                                              *
 * linked list functions                                         *
 *                                                               *
 * started 6/15/98, Karl Kosack (kosack@andrew.cmu.edu)          *
 *****************************************************************/
#ifndef __LIST_H__
#define __LIST_H__

#include <iostream.h>

template<class Type>
class ListItem{
	public:
	Type *node;
	ListItem *next;
};


template<class Type>
class List {

  public:
	List(void);
	~List(void);
	void  addItem(Type *);
	bool  isEmpty(void);
	void  beginIteration(void);
	Type* getItem(void);
	int  getLength(void);

	ListItem<Type> *current;
	ListItem<Type> *first;
  
};

template<class Type>
List<Type>::List(void) {
  first = NULL;
  current = first;
}

template<class Type>
List<Type>::~List(void){

	ListItem<Type> *temp;

	while (first != NULL && first->next != NULL){
		temp = first->next;		
		delete first;
		first = temp;
	}


	delete first;
	first = NULL;

}


template<class Type>
void
List<Type>::addItem(Type *item){
  ListItem<Type> *newitem;

  newitem = new ListItem<Type>;
  newitem->node = item;
  newitem->next = first;
  first = newitem;

}

template<class Type>
bool
List<Type>::isEmpty(void){
  if (first==NULL) return true;
  else return false;
}

template<class Type>
void
List<Type>::beginIteration(void){
  current = first;
}


template<class Type>
Type *
List<Type>::getItem(void){
  Type *temp;

  if(current == NULL) return NULL;
  temp = current->node;
  
  current = current->next;
  return temp;

}

template<class Type>
int
List<Type>::getLength(void){
	ListItem<Type>  *temp;
	int i=0;

	for(temp=first; temp != NULL; temp=temp->next){
		i++;
	}
	
	return i;

}

#endif //__LIST_H__