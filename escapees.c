#include "truth_assign.h"
#include <stdio.h>
#include <string.h>
#include "llist.h"
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"

#define ENG_RNG 0.001
#define COMP_INT 1667457891

/*
 * Looks to see if a gamma in a history escapes the detector module without
 * expending all of its energy. It does this by reading the history and looking
 * for a final photoelectric effect occurance in the volume. If no such
 * occurance occurs the number of scatters before leaving is printed to the
 * output file. If it is completely contained, a -1 is recorded instead.
 */

  
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
	// char origin[20];
	int count;

	int worked;

	worked = fscanf(source, "%u", &numb);
	worked = fscanf(source, "%lf", &energy);
	worked = fscanf(source, "%lf", &deposit);
	worked = fscanf(source, "%lf", &x);
	worked = fscanf(source, "%lf", &y);
	worked = fscanf(source, "%lf", &z);
	worked = fscanf(source, "%lf", &tof);
	worked = fscanf(source, "%i", &particle);
	// worked = fscanf(source, "%s", origin);
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
	new_event->deposited 	= deposit;
	new_event->location 	= three_vec(x,y,z);
	new_event->tof 			= tof;
	new_event->particle 	= particle;
	// strcpy(new_event->orgin, origin);
	new_event->id		= count;
	new_event->orgin[0]		= (char)0;


	// if (READ_DEBUG) {
	// 	print_event((void*)new_event);
	// }

	return new_event;
}

// frees event malloc, returns NULL. For fmap
void* delete_event(void* in) {
	if (in == NULL) {
		return NULL;
	}
	free(((event*)in)->location); // frees the allocated vector
	free(in);
	return NULL;
}

// just free with a return value of NULL for fmap
void* free_null(void* in) {
	if (in == NULL) {
		return NULL;
	}
	free(in);
	return NULL;
}

// prints an event, returns the event pointer
void* print_event(void* in) {
	event *val = (event*)in;
	printf("History: %u, Energy: %.1f, Particle type %i, ", val->number, val->energy, val->particle);
	printf("Identifier: %i.\n", val->id);
	return (void*)val;
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
		// history_num = 0;
		previous_event = f(source);
		if (previous_event == NULL) {
			// probably at EOF, in any case we need to be done
			return NULL;
		}
	}
	history_num = previous_event->number;
	while (history_num == previous_event->number) {
		history = add_to_bottom(history, previous_event);
		previous_event = f(source);
		if (previous_event == NULL) {
			// EOF reached
			return history;
		}
	}
	// make a new copy of the event in previous event for storing
	// means it will continue pointing right as otherwise it can
	// point to the history that got freed
	previous_event = duplicate_event(previous_event);
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
 * is_contained: Finds if the history has the gamma ending with no energy,
 * therefore the gamma is contained.
 */
int is_contained(llist* history) {
	while (history != NULL) {
		event* current_event = (event*)history->data;
		if (current_event->particle == 11) {
			// if ((current_event->orgin[0] == 'p') && (current_event->orgin[1] == 'h')) {
			return 1;
			// }
		}
		history = history->down;
	}
	return 0;
}

/*
 * total_scatters
 * counts the number of occurances of particles from the "Primary" origin in the
 * history. The count starts at -1 so that the entrance appearance does not
 * effect the count.
 */
int total_scatters(llist* history) {
	history = list_head(history);
	int total = -1;
	while (history != NULL) {
		event* current_event = (event*)history->data;
		if (current_event->id == 1) {
			total++;
		}
		history = history->down;
	}
	return total;
}

/*
 * exit_energy
 * returns the energy of the final primary generation gamma that appears in the
 * history. This should be the energy of the exit gamma if the gamma exited
 */
double exit_energy(llist* history) {
	history = list_head(history);
	event* most_recent_gamma = NULL;
	while (history != NULL) {
		event* current_event = (event*)history->data;
		if (current_event->orgin[0] == 'P') {
			if (current_event->particle == 22) {
				most_recent_gamma = current_event;
			}
		}
		history = history->down;
	}
	if (most_recent_gamma == NULL) {
		return -1;
	}
	return most_recent_gamma->energy;
}

int main(int argc, char **argv) {
	// quick reminder that argc is the number of arguments,
	// argv points to the strings of the arguments

	// we are looking for an input histories file, an output file name,
	// and an inner radius of the detector this may later be changed to
	// input histories, input data file, output file name
	if (argc != 3) {
		printf("Unable to run. Expected 2 arguments got %i.\n", argc - 1);
		printf("Expects an input volume history file");
		printf(" and the output file name.");
		return 1;
	}
	FILE* in_histories = fopen(argv[1], "r");
	if (in_histories == NULL) {
		printf("Unable to open history file\n");
		return 1;
	}

	FILE* out_containment = fopen(argv[2], "w");
	if (out_containment == NULL) {
		printf("Unable to open output file for writing\n");
		return 1;
	}


	llist *history = load_history(in_histories, read_line);
	// some extra vars for calculating percentage of in paitent scatters
	uint contained		= 0;
	uint hist_num		= 0;
	while (history != NULL) {
		
		if (is_contained(history)) {
			fprintf(out_containment, "%i, -1, %f\n", hist_num, -1.0);
			contained++;
		} else {
			fprintf(out_containment, "%i, ", hist_num);
			fprintf(out_containment, "%i, ", total_scatters(history));
			fprintf(out_containment, "%f\n", exit_energy(history));
		}


		// cleanup of the current tick
		fmap(history, delete_event);
		delete_list(history);
		history = load_history(in_histories, read_line);
		hist_num++;

	}

	printf("Contained:	%i\n", contained);
	printf("Total:		%i\n", hist_num);
	

	fclose(in_histories);
	fclose(out_containment);
	return 0;
}

