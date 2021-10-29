#include <stdio.h>
#include <stdlib.h>
#include "llist.h"

// adds a piece of data to the top of a list
llist* add_to_top(llist* list, void* data) {
	// does the list even exist yet?
	if (list == NULL) {
		// ok time to make a list
		list = (llist*)malloc(sizeof(llist));
		list->up = NULL;
		list->down = NULL;
		list->data = data;
		return list;
	}
	// ok the list exists, now get to top, then add an entry
	list = list_head(list);
	llist *new_top = (llist*)malloc(sizeof(llist));
	new_top->down = NULL;
	new_top->up = list;
	new_top->data = data;
	list->down = new_top;
	return list;
}

// returns the value of a given list component
void* get_value(llist* list, unsigned int id) {
	if (list == NULL) {
		return NULL;
	}
	list = list_head(list);
	if (id >= list_length(list)) {
		fprintf(stderr, "get_value: cannot get value of location past end of list");
		return NULL;
	}
	if (id == 0) {
		return list->data;
	}
	return get_value(list->down, id--);
}

// adds a piece of data to the end of a list
llist* add_to_bottom(llist* list, void* data) {
	// does the list even exist yet?
	if (list == NULL) {
		// ok time to make a list
		list = (llist*)malloc(sizeof(llist));
		list->down = NULL;
		list->up = NULL;
		list->data = data;
		return list;
	}

	// ok the list exists, now get to bottom, then add an entry
	list = list_tail(list);
	llist *new_tail = (llist*)malloc(sizeof(llist));
	new_tail->down = NULL;
	new_tail->up = list;
	new_tail->data = data;
	list->down = new_tail;
	return list_head(list);
}

int rec_delete_list(llist* list) {
	// takes a list, calls itself on the component below,
	// then deletes active place
	if (list == NULL) {
		return 1;
	}
	if (list->data != NULL) {
		// to avoid memory leak will send error message if data not NULL
		fprintf(stderr, "delete_list error: deleted list not NULL");
		rec_delete_list(list->down);
		free(list);
		return 0;
	}
	int success = rec_delete_list(list->down);
	free(list);
	return success;
}

// deletes an empty list, cannot clear the data
int delete_list(llist* list) {
	if (list == NULL) {
		return 1;
	}
	list = list_head(list);
	return rec_delete_list(list);
}

// returns the head of a list
llist* list_head(llist* list) {
	if (list == NULL) {
		return NULL;
	}
	if (list->up == NULL) {
		return list;
	}
	return list_head(list->up);
}

// returns the tail of a list
llist* list_tail(llist* list) {
	if (list == NULL) {
		return NULL;
	}
	if (list->down == NULL) {
		return list;
	}
	return list_tail(list->down);
}

// returns the length of a list
int list_length(llist* list) {
	if (list == NULL) {
		return 0;
	}
	list = list_head(list);
	//now at top of list
	unsigned int i = 1;
	while (list->down != NULL) {
		list = list->down;
		i++;
	}
	return i;
}


void fmap_rec(llist* list, void* (*f)(void*)) {
	if (list == NULL) {
		return;
	}
	list->data = f(list->data);
	// fmap(list->up, f);
	fmap_rec(list->down, f);
	return;
}

// applies the function f to every instance of data in the list
void fmap(llist* list, void* (*f)(void*)) {
	if (f == NULL) {
		return;
	}
	list = list_head(list);
	fmap_rec(list, f);
	return;
}

int rec_list_and(llist* list, int (*f)(void*)) {
	if (list == NULL) {
		return 1;
	}
	return f(list->data) && rec_list_and(list->down, f);
}

int list_and(llist* list, int (*f)(void*)) {
	if (f == NULL) {
		return 0;
	}
	list = list_head(list);
	return rec_list_and(list, f);
}

int rec_list_or(llist* list, int (*f)(void*)) {
	if (list == NULL) {
		return 0;
	}
	return f(list->data) || rec_list_or(list->down, f);
}

int list_or(llist* list, int (*f)(void*)) {
	if (f == NULL) {
		return 0;
	}
	list = list_head(list);
	return rec_list_or(list, f);
}