#include "truth_assign.h"
#include <stdio.h>
#include <string.h>
#include "llist.h"
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"

#define ENG_RNG 1.
#define COMP_INT 1667457891


  
// reads a line from the source file as an event
event* read_line(FILE* source) {

	// begin to fill the event struct with information from .phsp
	uint numb;
	double energy;
	double deposit;
	double x;
	double y;
	double z;
	double tof;
	int particle;
	char origin[20];
	int count;

	int worked;

	fscanf(source, "%u", &numb);
	fscanf(source, "%lf", &energy);
	fscanf(source, "%lf", &deposit);
	fscanf(source, "%lf", &x);
	fscanf(source, "%lf", &y);
	fscanf(source, "%lf", &z);
	fscanf(source, "%lf", &tof);
	fscanf(source, "%i", &particle);
	fscanf(source, "%s", origin);
	worked = fscanf(source, "%i", &count);

	if (worked == EOF) {
		return NULL;
	}
	
	// make a new event to be passed out
	event* new_event = (event*)malloc(sizeof(event));
	if (new_event == NULL) {
		return NULL;
	}
	new_event->number 		= numb;
	new_event->energy 		= energy;
	new_event->depoisted 	= deposit;
	new_event->location 	= three_vec(x,y,z);
	new_event->tof 			= tof;
	new_event->particle 	= particle;
	strcpy(new_event->orgin, origin);
	new_event->id		= count;

	return new_event;
}

// frees event malloc, returns NULL. For fmap
void* delete_event(void* in) {
	event *val = (event*)in;
	free(val->location);
	free(val);
	return NULL;
}

// prints an event, returns the event pointer
void* print_event(void* in) {
	event *val = (event*)in;
	printf("History: %u, Energy: %.1f, Particle type %i, ", val->number, val->energy, val->particle);
	printf("Identifier: %i.\n", val->id);
	return (void*)val;
}

void* print_event_to_out (void* in) {
	event* val = (event*)in;
	fprintf(stdout, "%12u ", val->number);
	fprintf(stdout, "%12f ", val->energy);
	fprintf(stdout, "%12f ", val->depoisted);
	fprintf(stdout, "%12f ", val->location->x);
	fprintf(stdout, "%12f ", val->location->y);
	fprintf(stdout, "%12f ", val->location->z);
	fprintf(stdout, "%12f ", val->tof);
	fprintf(stdout, "%12i ", val->particle);
	fprintf(stdout, "%22s ", val->orgin);
	fprintf(stdout, "%12i ", val->id);
	fprintf(stdout, "\n");
}

// loads all of the events in a history.
// also ends up getting the first event of the next history
// to load it uses a function that take the file type and
// outputs an event*
llist* load_history(FILE* source, event* (*f)(FILE*)) {

	static event* previous_event = NULL;
	uint history_num;
	llist* history = NULL;
	if (previous_event == NULL) {
		history_num = 0;
		previous_event = f(source);
		if (previous_event == NULL) {
			// probably at EOF, in any case we need to be done
			return NULL;
		}
	} else {
		history_num = previous_event->number;
	}
	while (history_num == previous_event->number) {
		history = add_to_bottom(history, previous_event);
		previous_event = f(source);
		if (previous_event == NULL) {
			// EOF reached
			return history;
		}
	}
	return history;
}

// loads all of the events in a history.
// also ends up getting the first event of the next history
// to load it uses a function that take the file type and
// outputs an event*
llist* load_historyb(FILE* source, event* (*f)(FILE*)) {

	static event* previous_event = NULL;
	uint history_num;
	llist* history = NULL;
	if (previous_event == NULL) {
		history_num = 0;
		previous_event = f(source);
		if (previous_event == NULL) {
			// probably at EOF, in any case we need to be done
			return NULL;
		}
	} else {
		history_num = previous_event->number;
	}
	while (history_num == previous_event->number) {
		history = add_to_bottom(history, previous_event);
		previous_event = f(source);
		if (previous_event == NULL) {
			// EOF reached
			return history;
		}
	}
	return history;
}

/* 
 * double_equality:
 * determines if two doubles are within range of each other. if they
 * are, then it returns 1, otherwise it returns zero
 */
int double_equality(double a, double b, double range) {
	double diff = a - b;
	if ((diff >= -range) && (diff <= range)) {
		return 1;
	}
	return 0;
}

/* 
 * Takes an event* and makes a second event for placement elsewhere. Uses malloc
 * to make a new event, then copies the information over
 */
event* copy_event(event* source) {
	if (source == NULL) {
		return NULL;
	}
	event* new_event = (event*)malloc(sizeof(event));
	if (new_event == NULL) {
		fprintf(stderr, "copy_event: malloc failed");
		return NULL;
	}
	new_event->number 		= source->number;
	new_event->energy 		= source->energy;
	new_event->depoisted 	= source->depoisted;
	new_event->location 	= three_vec(source->location->x,
										source->location->y,
										source->location->z);
	new_event->tof			= source->tof;
	new_event->particle		= source->particle;
	strcpy(new_event->orgin, source->orgin);
	new_event->id 			= source->id;
	fprintf(stderr, "hello world\n");
	return new_event;
}

int main(int argc, char **argv) {
	// quick reminder that argc is the number of arguments,
	// argv points to the strings of the arguments

	// we are looking for an input histories file, an output file name,
	// and an inner radius of the detector this may later be changed to
	// input histories, input data file, output file name
	if (argc != 3) {
		printf("Unable to run. Expected 2 arguments got %i.\n", argc - 1);
		printf("Expects an input full history file and input virtual volume  ");
		printf("history file\n");
		// printf(" inner radius, detector half height and energy cutoff value.\nDetector sizes");
		// printf(" are in the same units as the input history file's distances.\n");
		return 1;
	}
	FILE* in_histories = fopen(argv[1], "r");
	if (in_histories == NULL) {
		printf("Unable to open history file\n");
		return 1;
	}
	FILE* virt_file = fopen(argv[2], "r");
	if (virt_file == NULL) {
		printf("Unable to open detector history file\n");
		return 1;
	}


	llist *history = load_history(in_histories, read_line);
	llist *vert_history = load_historyb(virt_file, read_line);
	// some extra vars for calculating percentage of in paitent scatters
	uint hist_num = 0;
	while (history != NULL) {
		// if the virtual history have more than one entry, check that the last
		// entry of the virtual history has about the same energy as the last
		// entry of the real volume list. If it does, then append the virtual
		// history event to the real history and print the new real history using
		// the print print_event_to_out function.

		if (list_length (vert_history) > 1) {
			// check if the last entry of the real and virtual history have
			// about the same energy (as expected for a common gamma)
			double hist_energy = ((event*)list_tail(history))->energy;
			event* virt_tail = (event*)(list_tail(vert_history)->data);
			if (double_equality(hist_energy, virt_tail->energy, ENG_RNG)) {
				add_to_bottom(history, copy_event(virt_tail));
			}
		}
		fmap(history, print_event_to_out);



		fmap(history, delete_event);
		delete_list(history);
		history = load_history(in_histories, read_line);
		hist_num++;
		// load a new history for the second list if the second list is out of date
		if (history != NULL) {	
			// check that there is still information left in the main history
			if (((event*)history->data)->number > ((event*)vert_history->data)->number) {
				// vert_history is out of date
				fmap(vert_history, delete_event);
				delete_list(vert_history);
				vert_history = load_historyb(virt_file, read_line);
			}
		}

	}


	fclose(in_histories);
	fclose(virt_file);
	return 0;
}

