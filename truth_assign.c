#include "truth_assign.h"
#include <stdio.h>
#include <string.h>
#include "llist.h"
#include <math.h>
#include <stdlib.h>

#define ENG_RNG 0.001
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
	new_event->location[0] 	= x;
	new_event->location[1] 	= y;
	new_event->location[2] 	= z;
	new_event->tof 			= tof;
	new_event->particle 	= particle;
	strcpy(new_event->orgin, origin);
	new_event->id		= count;

	return new_event;
}

// frees event malloc, returns NULL. For fmap
void* delete_event(void* in) {
	event *val = (event*)in;
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
 * Checks to see if the current event was a gamma event (scattering)
 * closer to the center than the given distance. If it was it returns
 * 1 (true). If it was not inside it returns a
 * zero (false). Designed to be fast (so that it can be run on every
 * event).
 */
uint inside_radius(double distance, double height, event* event) {
	if (event->particle != 22) {
		// not a gamma, can't be in paitent scatter
		// printf("not a gamma\n");
		return 0;
	}
	double radius2 = event->location[0] * event->location[0] + event->location[1] * event->location[1];
	// printf("particle %i, distance^2 %d", event->particle, event->count);
	if (sqrt(radius2) < distance) {
		if (fabs(event->location[2]) < height){
			return 1;
		}
		return 0;
	}
	return 0;
}

int rec_in_paitent(llist* list, double radius, double height) {
	if (list == NULL) {
		return 0;
	}
	if (!inside_radius(radius, height, (event*)(list->data))) {
		// not inside in paitent radius, can't be in paitent scatter
		// printf("rec_in_paitent: not inside\n");
		return rec_in_paitent(list->down, radius, height);
	}
	// we now have a mention of a gamma inside the in paitent
	// radius, now to check if it lost energy from it's previous
	// mention
	if ((list->up != NULL) && (((event*)(list->data))->id == ((event*)(list->up->data))->id)) {
		// there is a previous mention of the same particle
		// printf("Particle in loc w/ predicessor\n");
		if (((event*)(list->data))->energy < ((event*)(list->up->data))->energy) {
			// ok now the previous interaction of the same particle had a larger
			// energy by more than the chosen energy range. We have an in paitent scatter!
			// we can stop looking now.
			event* ev = (event*)list->data;
			// printf("Found in paitent scattering, history %i, particle %i, ", ev->number, ev->count);
			// printf("pre-scatter energy %f, post-scatter energy %f\n", ((event*)(list->up->data))->energy, ev->energy);
			return 1;
		// } else {
		// 	event* ev = (event*)list->data;
		// 	printf("pre-scatter energy %f, post-scatter energy %f\n", ((event*)(list->up->data))->energy, ev->energy);
		}
	}
	// well better luck next time
	// printf("rec_in_paitent: other problem\n");
	return rec_in_paitent(list->down, radius, height);
}

/* 
 * Checks to see if there has been an in paitent scattering event in
 * the given history. An event is defined as a gamma showing up in the
 * list within the given radius of the z axis, with the previous instance
 * of the gamma having an energy greater by at least ENG_RNG. It finds
 * this information by calling rec_in_paitent, which checks over the
 * list for a gamma that fits the above conditions
 */
int in_patient(llist* list, double radius, double height) {
	list = list_head(list);
	return rec_in_paitent(list, radius, height);
}



int main(int argc, char **argv) {
	// quick reminder that argc is the number of arguments,
	// argv points to the strings of the arguments

	// we are looking for an input histories file, an output file name,
	// and an inner radius of the detector this may later be changed to
	// input histories, input data file, output file name
	if (argc != 5) {
		printf("Unable to run. Expected 4 arguments got %i.\n", argc - 1);
		printf("Expects an input history file, output file name, detector");
		printf(" inner radius and detector half height.\nDetector sizes");
		printf(" are in the same units as the input history file's distances.\n");
		return 1;
	}
	FILE* in_histories = fopen(argv[1], "r");
	if (in_histories == NULL) {
		printf("Unable to open history file\n");
		return 1;
	}
	FILE* out_in_patient = fopen(argv[2], "w");
	if (out_in_patient == NULL) {
		printf("Unable to open output file for writing\n");
		return 1;
	}
	double in_paitent_distance = strtod(argv[3], NULL);
	double detector_height = strtod(argv[4], NULL);
	if ((in_paitent_distance <= 0) || (detector_height <= 0)) {
		printf("Size dimensions must be greater than zero\n");
		return 1;
	}

	// current loop version, just for testing
	llist *history = load_history(in_histories, read_line);
	// some extra vars for calculating percentage of in paitent scatters
	uint numerator = 0;
	uint denominator = 0;
	while (history != NULL) {
		// fmap(history, print_event);
		if (in_patient(history, in_paitent_distance, detector_height)) {
			fprintf(out_in_patient, "%i, true\n", ((event*)(history->data))->number);
			numerator++;
			denominator++;
		} else {
			fprintf(out_in_patient, "%i, false\n", ((event*)(history->data))->number);
			denominator++;
		}
		fmap(history, delete_event);
		delete_list(history);
		history = load_history(in_histories, read_line);
	}
	if (denominator != 0) {
		printf("Percentage of in paitent scatters: ");
		float percent = (float)numerator / (float)denominator;
		printf("%f\n", percent);
	}


	fclose(in_histories);
	fclose(out_in_patient);
	return 0;
}

