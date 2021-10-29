#ifndef llist_h
#define llist_h



typedef struct llist_ {
	struct llist_* up;
	void* data;
	struct llist_* down;
} llist;

// adds a piece of data to the top of a list
llist* add_to_top(llist* list, void* data);

// returns the value of a given list component
void* get_value(llist* list, unsigned int id);

// adds a piece of data to the end of a list
llist* add_to_bottom(llist* list, void* data);

// deletes an empty list, cannot clear the data
int delete_list(llist* list);

// returns the head of a list
llist* list_head(llist* list);

// returns the tail of a list
llist* list_tail(llist* list);

// returns the length of a list
int list_length(llist* list);

// applies the function f to every instance of data in the list
void fmap(llist* list, void* (*f)(void*));

// ands the result of f across the whole list
int list_and(llist* list, int (*f)(void*));

// ors the result of f across the whole list
int list_or(llist* list, int (*f)(void*));

#endif