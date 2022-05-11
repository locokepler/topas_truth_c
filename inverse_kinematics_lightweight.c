#include "truth_assign.h"
#include <stdio.h>
#include <string.h>
#include "llist.h"
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"
#include <pthread.h>

#define MAX_THREAD_CALLS 8

#define ENG_RNG 0.001f
#define COMP_INT 1667457891
#define ELECTRON_MASS 510.999f
#define SPD_LGHT 29.98f
#define LARGEST 10

#define DEP_UNCERT 5.f
#define SPC_UNCERT 0.5f


#define READ_DEBUG 0
#define GENERAL_DEBUG 0

pthread_t tid[MAX_THREAD_CALLS];
pthread_mutex_t file_lock;

  
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

	worked = fscanf(source, "%u", &numb);
	worked = fscanf(source, "%lf", &energy);
	worked = fscanf(source, "%lf", &deposit);
	worked = fscanf(source, "%lf", &x);
	worked = fscanf(source, "%lf", &y);
	worked = fscanf(source, "%lf", &z);
	worked = fscanf(source, "%lf", &tof);
	worked = fscanf(source, "%i", &particle);
	worked = fscanf(source, "%s", origin);
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

	// if (READ_DEBUG) {
	// 	print_event((void*)new_event);
	// }

	return new_event;
}

// creates a new scatter structure and fills it
scatter* new_scatter_old(vec3d* vector, double deposited, double time) {
	scatter* new = (scatter*)malloc(sizeof(scatter));
	new->deposit = deposited;
	new->loc = vector;
	new->dir = NULL;
	new->time = time;
	new->eng_uncert = -1;
	new->space_uncert = -1;
	new->time_uncert = -1;
	return new;
}

scatter* new_scatter(vec3d* vector, vec3d* dir, double deposit, double time, double eng_uncert, double space_uncert, double time_uncert) {
	scatter* new = (scatter*)malloc(sizeof(scatter));
	new->deposit = deposit;
	new->eng_uncert = eng_uncert;
	new->space_uncert = space_uncert;
	new->loc = vector;
	new->dir = dir;
	new->time = time;
	new->time_uncert = time_uncert;
	new->truth = NULL;
	return new;
}

scatter* copy_scatter(scatter* a) {
	return new_scatter(vec_copy(a->loc),vec_copy(a->dir), a->deposit, a->time, a->eng_uncert, a->space_uncert, a->time_uncert);
}

void* delete_scatter(void* in) {
	if (in == NULL) {
		return NULL;
	}
	free(((scatter*)in)->loc);
	free(((scatter*)in)->dir);
	if (((scatter*)in)->truth != NULL) {
		free(((scatter*)in)->truth);
	}
	free(in);
	return NULL;
}

event* duplicate_event(event* source) {
	event* new_event = (event*)malloc(sizeof(event));
	new_event->number		= source->number;
	new_event->energy		= source->energy;
	new_event->depoisted	= source->depoisted;
	new_event->location		= vec_copy(source->location);
	new_event->tof			= source->tof;
	new_event->particle		= source->particle;
	strcpy(new_event->orgin, source->orgin);
	new_event->id			= source->id;
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

void* free_lor(void* in) {
	if (in == NULL) {
		return NULL;
	}
	lor* clear = (lor*)in;
	free(clear->center);
	free(clear->dir);
	free(clear);
	return NULL;
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
			if (GENERAL_DEBUG) {
				printf("NULL event reached, ending\n");
			}
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
			if (GENERAL_DEBUG) {
				printf("NULL event reached, ending.\n");
				printf("History number %i\n", history_num);
			}
			return history;
		}
	}
	// make a new copy of the event in previous event for storing
	// means it will continue pointing right as otherwise it can
	// point to the history that got freed
	previous_event = duplicate_event(previous_event);
	if (GENERAL_DEBUG) {
		printf("History number %i\n", history_num);
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

// now for a quick array sorting implemetation

int scatter_dep_compare(void* va, void* vb) {
	if (va == NULL || vb == NULL) {
		return 0;
	}
	scatter* a = (scatter*)va;
	scatter* b = (scatter*)vb;
	if (a->deposit > b->deposit) {
		return 1;
	} else if (a->deposit == b->deposit) {
		return 0;
	}
	return -1;
}

int partition(void** array, int low, int high, int (*f)(void*, void*)) {
	if (array == NULL) {
		return 0;
	}
	void* pivot = array[high];
	// the pivot 
	
	int i = low - 1;
	// the current best location for the pivot
	for (int j = low; j < high - 1; j++) {
		if (f(array[j], pivot) < 0) {
			i++;
			void** holding = array[j];
			array[j] = array[i];
			array[i] = holding;
		}
	}
	void** holding = array[i + 1];
	array[i + 1] = array[high];
	array[high] = holding;
	return i + 1;
}

void scatter_quicksort(scatter** arr, int low, int high) {
	if (arr == NULL) {
		return;
	}
	if (low < high) {
		// pivot
		int pivot = partition((void**)arr, low, high, scatter_dep_compare);

		scatter_quicksort(arr, low, pivot -1);
		scatter_quicksort(arr, pivot + 1, high);
	}
}
// finish array sorting implementation

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
	double radius2 = event->location->x * event->location->x 
						+ event->location->y * event->location->y;
	// printf("particle %i, distance^2 %d", event->particle, event->count);
	if (sqrt(radius2) < distance) {
		if (fabs(event->location->z) < height){
			return 1;
		}
		return 0;
	}
	return 0;
}

int rec_in_paitent(llist* list, double radius, double height, int id) {
	if (list == NULL) {
		return 0;
	}
	if (((event*)list->data)->id != id) {
		return rec_in_paitent(list->down, radius, height, id);
		// not our particle
	}
	if (!inside_radius(radius, height, (event*)(list->data))) {
		// not inside in paitent radius, can't be in paitent scatter
		// printf("rec_in_paitent: not inside\n");
		return rec_in_paitent(list->down, radius, height, id);
	}
	// we now have a mention of the gamma inside the in paitent
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
			// even though we have an in patient scatter we want to check we don't
			// have a second one or more from this gamma
			if (ev->energy < energy_cutoff)
				return 100; 
				// returns a large number of scatters so that the value
				// can be spotted and rejected later in the code.
			return 1 + rec_in_paitent(list->down, radius, height, id);
			// will return >1 if an additional in patient scatter occcurs
		}
	}
	// well better luck next time
	// printf("rec_in_paitent: other problem\n");
	return rec_in_paitent(list->down, radius, height, id);
}

/* 
 * Checks to see if there has been an in paitent scattering event in
 * the given history. An event is defined as a gamma showing up in the
 * list within the given radius of the z axis, with the previous instance
 * of the gamma having an energy greater by at least ENG_RNG. It finds
 * this information by calling rec_in_paitent, which checks over the
 * list for a gamma that fits the above conditions
 */
int in_patient(llist* list, double radius, double height, int id) {
	if (list == NULL)
		fprintf(stderr, "in_patient passed NULL pointer\n");
	list = list_head(list);
	return rec_in_paitent(list, radius, height, id);
}

/*
 * find_annihilation_point
 * Takes a history and finds where the annihilation occured. If no
 * annihilation occured it returns NULL. If one did occur it returns an
 * array of doubles of length 3. This has the x,y,z position of the final
 * location of the positron
 */
vec3d* find_annihilation_point(llist *history) {
	if (history == NULL) {
		return NULL;
	}
	history = list_head(history);
	// looks for the first event that refers to a positron
	while (((event*)(history->data))->particle != -11) {
		history = history->down;
		if (history == NULL) {
			return NULL;
		}
	}
	// now look for the last event that refers to a positron
	while ((history->down != NULL) && (((event*)(history->down->data))->particle == -11)) {
		history = history->down;
	}
	// now at the last event, return the location as a 3-vector
	return ((event*)(history->data))->location;
}

/*
 * line_to_dot_dist:
 * finds the distance from a line defined by the start and end points to a
 * given point. The distance is how far the minimum distance is.
 */
double line_to_dot_dist(vec3d* start, vec3d* end, vec3d* point) {
	vec3d* num_first_term = vec_sub(start, end);
	vec3d* num_sec_term = vec_sub(start, point);
	// vec3* denom_vec = vec_sub(end, start);
	double numerator = vec_mag(vec_cross(num_first_term, num_sec_term));
	double denomenator = vec_mag(num_first_term);
	free(num_first_term);
	free(num_sec_term);
	// free(denom_vec);
	return numerator / denomenator;
}

// does what it says on the tin. Adds a and b in quadrature
double add_quadrature(double a, double b) {
	return sqrt((a * a) + (b * b));
}

// using the partial derviative in b and theta computes the uncertanty in the
// energy of the gamma going to b
double expected_uncert_b(double b, double theta, double uncert_b, double uncert_theta) {
	// if given bullshit uncertanties, say it
	if ((uncert_b < 0) || (uncert_theta < 0)) {
		return -1.;
	}

	// first partial derivative in theta
	double numerator = 511. * b * sin(theta);
	double one_less_theta = 1. - cos(theta);
	double root = b * b + ((4. * 511. * b) / one_less_theta);
	double denom = one_less_theta * one_less_theta * sqrt(root);
	double theta_der = numerator / denom;
	uncert_theta = uncert_theta * theta_der;

	// now partial derivative in b
	numerator = b + ((2. * 511.) / one_less_theta);
	root = b * b + ((4. * 511. * b) / one_less_theta);
	double first_term = numerator / sqrt(root);
	double b_der = .5 * (first_term = 1.);
	uncert_b = uncert_b * b_der;

	return add_quadrature(uncert_theta, uncert_b);
}

/* 
 * expected_energy_b
 * takes 3 scatters and solves the kinematics assuming that the gamma goes from
 * a->b->c, solving for the energy of the gamma from a->b. To do this it
 * calculates the angle <ABC, then uses the deposited energy at B and this angle
 * to determine the energy of the incoming gamma. 
 */
 double expected_energy_b(scatter* a, scatter* b, scatter* c, double* uncert) {
	 if ((a == NULL) || (b == NULL) || (c == NULL)) {
		 return 1.;
	 }

	// first calculate the angle at b

	// get the vector from b->a
	vec3d* ab = vec_sub(b->loc, a->loc);
	// get the vector from b->c
	vec3d* bc = vec_sub(c->loc, b->loc);
	// calculate the angle itself
	double theta = vec_angle(ab, bc);

	// calculate an estimate of uncertanty in theta
	// done by sin(theta) = theta for small angles
	double ab_uncert = add_quadrature(a->space_uncert, b->space_uncert);
	double bc_uncert = add_quadrature(c->space_uncert, b->space_uncert);
	ab_uncert = ab_uncert / vec_mag(ab);
	bc_uncert = bc_uncert / vec_mag(bc);
	double theta_uncert = add_quadrature(ab_uncert, bc_uncert);

	if (1 - cos(theta) == 0.) {
		// if this were the case things would break in division
		return -1.; // we can't tell what energy the gamma had, so return easy
		// to spot garbage
	}

	// free the used vectors
	free(ab);
	free(bc);


	// calculate the incoming gamma energy. 
	// The calculation is in the form of the quadradic formula. The first term
	// calculated is -4ac
	double last_term = 4 * ((b->deposit * 511.) / (1. - cos(theta)));
	double determinator = (b->deposit * b->deposit) + last_term;
	if (determinator < 0.) {
		// we don't want to take the square root of a negative number,
		// also it appears this scatter was non-physical. Honestly how did you
		// get here? you added two (theoretically) positive numbers and got a
		// negative one
		return -1.;
	}
	double gamma_to_b_e = (b->deposit + sqrt(determinator)) / 2.; // the energy of the
	// gamma from a to b

	// now to calculate the uncertanty
	if (uncert != NULL) {
		double delta_e = expected_uncert_b(b->deposit, theta, b->eng_uncert, theta_uncert);
		*uncert = delta_e;
	}

	return gamma_to_b_e;
 }

/* 
 * expected_energy_a
 * Takes three locations given as vec3 pointers. Energies are given as doubles.
 * The energies are assumed to be in keV.
 * The output is the estimated energy of the gamma before interacting at a.
 * This function is designed for use in hypothesis testing.
 * Method: First it computes the ange <ABC. Then, using the energy deposited by
 * the Compton scattering at B and the angle caused it calcuates the incoming
 * gamma's energy. The energy deposited at A is then added to the gamma energy
 * coming to B giving the energy of the incident gamma to A.
 */
double expected_energy_a(scatter* a, scatter* b, scatter* c) {
	double gamma_to_b_e = expected_energy_b(a, b, c, NULL);

	if (gamma_to_b_e < 0) {
		return gamma_to_b_e;
		// covers error codes, all with negative energy returns. This includes
		// if a is null, as then expected_energy_b cannot run, returns -1.
	}

	// now just add the energy of the gamma a->b and the energy deposited at a
	return gamma_to_b_e + a->deposit;


}

/*
 * test_expected_energy
 * a test suite for the expected energy function to check that it is producing
 * the expected results.
 * Only has to test expected_energy_a, as expected_energy_b is used as the main
 * component of a.
 */
int test_expected_energy() {
	// test is a set of scatters with a known starting energy
	// input is a 511. keV gamma, scatters with 10. keV deposit at A
	// then scatters with an angle of 120 degrees at B, deposit 203.1654 keV
	// then goes to C.
	int pass = 1;

	vec3d* point_a = three_vec(0., 0., 0.);
	vec3d* point_b = three_vec(0., 3., 0);
	vec3d* point_c = three_vec(0., 2., 1.73205);
	double deposit_a = 127.405;
	double deposit_b = 203.1654;
	scatter* scatter_a = new_scatter_old(point_a, deposit_a, -1);
	scatter* scatter_b = new_scatter_old(point_b, deposit_b, -1);
	scatter* scatter_c = new_scatter_old(point_c, 10., -1);
	// run the function
	double result = expected_energy_a(scatter_a, scatter_b, scatter_c);
	delete_scatter(scatter_a);
	delete_scatter(scatter_b);
	delete_scatter(scatter_c);
	if (double_equality(result, 511., 0.1)) {
		printf("test_expected_energy passed test 1\n");
	} else {
		printf("test_expected_energy FAILED test 1:\n");
		printf("expected energy %f, got energy %f.\n", 511., result);
		pass = 0;
	}


	if (pass) {
		printf("test_expected_energy: all test passed\n");
		return 1;
	}
	printf("test_expected_energy: FAILED ONE OR MORE TESTS\n");
	return 0;
}

/* 
 * scattering_iterator:
 * takes a history of scatters from a detector (position and energy information)
 * and determines the event it expects to be the true first scatter. It then
 * returns this scatter.
 * history->data must be of type scatter*
 * Method: Calls expected_energy_a on every possible orientation of the list of
 * scatters. If a scatter orientation has a closer energy to 511. keV then that
 * orientation is kept. If not it is ignored. Then the next orientation it run.
 */
scatter* scattering_iterator(llist* scatter_list, double energy_percent) {
	if (scatter_list == NULL) {
		fprintf(stderr, "scattering_iterator: no history given\n");
		return NULL;
	}

	// the current best find energy
	double best_find = 0;
	// the current best find scatter
	scatter* best_scatter = NULL;
	// length of the scatter list
	int list_len = list_length(scatter_list);

	if (list_len == 1) {
		// only 1 scatter exists, return it
		return NULL;
		// return (scatter*)history->data;
	}

	if (list_len == 2) {
		// only 2 scatters exit. return the brighter one
		return NULL;
		// if (((scatter*)history->data)->deposit > ((scatter*)history->down->data)->deposit) {
		// 	return (scatter*)history->data;
		// } else {
		// 	return (scatter*)history->down->data;
		// }
	}

	// time to turn the scattering list into an array for fast access
	scatter** locations = (scatter**)malloc(list_len * sizeof(scatter*));
	for (int i = 0; i < list_len; i++)
	{
		locations[i] = (scatter*)scatter_list->data;
		scatter_list = scatter_list->down;
	}
	


	double hypoth;
	// double second_best = 0;
	// now to iterate over all of the loop possibilites
	for (int first = 0; first < list_len; first++)
	{
		for (int second = 0; second < list_len; second++)
		{
			for (int third = 0; third < list_len; third++)
			{
				if ((first != second) && (first != third) && (second != third)) {
					hypoth = expected_energy_a(locations[first],
											locations[second],
											locations[third]);
					// check if the hypothesis is better than previous
					if (fabs(hypoth - ELECTRON_MASS) < fabs(best_find - ELECTRON_MASS)) {
						// the current hypthesis is better than the previous best
						if (GENERAL_DEBUG) {
							printf("Scattering Iterator: new best scatter found:\n");
							printf("%f keV at ", hypoth);
							vec_print(locations[first]->loc, stdout);
							printf("\n");
						}

						// second_best = best_find;
						best_find = hypoth;
						best_scatter = locations[first];
					}
				}
			}	
		}
	}

	// if (first_scat_hypot == 0) {
	// 	first_scat_hypot = fabs(best_find - second_best);
	// } else {
	// 	second_scat_hypot = fabs(best_find - second_best);
	// }

	free(locations);
	// done iterating, now have the best scatter found in the list
	if (fabs(best_find - ELECTRON_MASS) < (energy_percent * ELECTRON_MASS)) {
		// the result was within the energy cut
		return best_scatter;
	}
	// no result found within the energy cutoff
	return NULL;
}

/*
 * takes two scattering histories. Searches for the best guess as to the first
 * scatter in the first history. To do this it uses the locations in the second
 * history as endpoints of the LOR, then proceeds to do a 3 point approximation
 * of the kinematics of the Compton scattering. It chooses the scatter with the
 * smallest deviation from the expected energy (by probability of deviation).
 * This is then combined with the probability of that scatter being the first
 * using the brightest scatter ordering. If multiple good canidates are found,
 * the iteration continues further down the chain looking for cases with no
 * likely solutions.
 * 
 * inputs:
 * history1: a list of scatters with location, deposited energy, uncertanty in
 * position, and uncertanty in deposited energy. This is the locations where
 * the first scatter is being seached for. Only runs if this has 2 or more
 * scatters in it
 * 
 * history2: same type as history1, just the locations where the other side of
 * the LOR should be. Only runs if this has 1 or more scatters in it.
 * 
 * energy_percent: what the percentage cut is for acceptable scattering
 */
scatter* multi_gamma_iterator(llist* history1, llist* history2, double energy_percent) {
	if ((history1 == NULL) || (history2 == NULL)) {
		fprintf(stderr, "multi_gamma_iterator: empty history given");
		return NULL;
	}

	// first is setup. Best scatter found, behavior for undersize histories,
	// copy histories into arrays for fast access.

	double best_find = 0;

	scatter* best_scatter = NULL;

	int len_hist1 = list_length(history1);
	int len_hist2 = list_length(history2);

	if (len_hist1 < 2) {
		// history1 is too short to run the iteration process.
		return NULL;
	}
	if (len_hist2 < 1) {
		// idk if it is possible to even get here, but we can't continue if we do
		return NULL;
	}

	// time to turn the scattering list into an array for fast access
	scatter** scatters1 = (scatter**)malloc(len_hist1 * sizeof(scatter*));
	for (int i = 0; i < len_hist1; i++)	{
		scatters1[i] = (scatter*)history1->data;
		history1 = history1->down;
	}

	scatter** scatters2 = (scatter**)malloc(len_hist2 * sizeof(scatter*));
	for (int i = 0; i < len_hist2; i++)	{
		scatters2[i] = (scatter*)history2->data;
		history2 = history2->down;
	}

	// a couple of variables for debug information
	double hypoth;
	double second_best = 0;


	// time to iterate over all of the possible combinations. The avaliable
	// configurations are:
	// for all i in len_hist2 and all of j != k in len_hist1
	for (int i = 0; i < len_hist2; i++) {
		for (int j = 0; j < len_hist1; j++) {
			for (int k = 0; k < len_hist1; k++) {
				// can only do 3 point checks with j and k not being the same
				if (j != k) {
					hypoth = expected_energy_b(scatters2[i], scatters1[j], scatters1[k], NULL);
					// check if the hypothesis is better than previous
					if (fabs(hypoth - ELECTRON_MASS) < fabs(best_find - ELECTRON_MASS)) {
						// the current hypthesis is better than the previous best
						if (GENERAL_DEBUG) {
							printf("2 hist Scattering Iterator: new best scatter found:\n");
							printf("%f keV at ", hypoth);
							vec_print(scatters1[j]->loc, stdout);
							printf("\n");
						}
						second_best = best_find;
						best_find = hypoth;
						best_scatter = scatters1[j];
					}
				}
			}
		}
	}

	if (first_scat_hypot == 0) {
		first_scat_hypot = fabs(best_find - second_best);
	} else {
		second_scat_hypot = fabs(best_find - second_best);
	}

	free(scatters1);
	free(scatters2);
	// done iterating, now have the best scatter found in the list
	if (fabs(best_find - ELECTRON_MASS) < (energy_percent * ELECTRON_MASS)) {
		// the result was within the energy cut
		if (GENERAL_DEBUG) {
			printf("multi_gamma_iterator: best scatter solution found at %f kev at:\n", best_find);
			vec_print(best_scatter->loc, stdout);
			printf("\n");
		}
		return best_scatter;
	}
	// no result found within the energy cutoff
	return NULL;
}

double scatter_dir_dot(scatter* first, scatter* second) {
	if (second->dir == NULL) {
		return -1;
	}
	vec3d* a_less_b = vec_sub(first->loc, second->loc);
	double dot = vec_dot(a_less_b, second->dir);
	free(a_less_b);
	return dot;
}

/*
 * takes two scattering histories. Searches for the best guess as to the first
 * scatter in the first history. To do this it uses the locations in the second
 * history as endpoints of the LOR, then proceeds to do a 3 point approximation
 * of the kinematics of the Compton scattering. It chooses the scatter with the
 * smallest deviation from the expected energy (by probability of deviation).
 * Any results that would requrire the electron to point in a non-physical
 * direction are rejected.
 * 
 * inputs:
 * history1: a list of scatters with location, deposited energy, uncertanty in
 * position, and uncertanty in deposited energy. This is the locations where
 * the first scatter is being seached for. Only runs if this has 2 or more
 * scatters in it
 * 
 * history2: same type as history1, just the locations where the other side of
 * the LOR should be. Only runs if this has 1 or more scatters in it.
 * 
 * energy_percent: what the percentage cut is for acceptable scattering
 */
scatter* multi_gamma_ele_iterator(llist* history1, llist* history2, double energy_percent) {
	if ((history1 == NULL) || (history2 == NULL)) {
		fprintf(stderr, "multi_gamma_ele_iterator: empty history given");
		return NULL;
	}

	// first is setup. Best scatter found, behavior for undersize histories,
	// copy histories into arrays for fast access.

	double best_find = 0;

	scatter* best_scatter = NULL;

	int len_hist1 = list_length(history1);
	int len_hist2 = list_length(history2);

	if (len_hist1 < 2) {
		// history1 is too short to run the iteration process.
		return NULL;
	}
	if (len_hist2 < 1) {
		// idk if it is possible to even get here, but we can't continue if we do
		return NULL;
	}

	// time to turn the scattering list into an array for fast access
	scatter** scatters1 = (scatter**)malloc(len_hist1 * sizeof(scatter*));
	for (int i = 0; i < len_hist1; i++)	{
		scatters1[i] = (scatter*)history1->data;
		history1 = history1->down;
	}

	scatter** scatters2 = (scatter**)malloc(len_hist2 * sizeof(scatter*));
	for (int i = 0; i < len_hist2; i++)	{
		scatters2[i] = (scatter*)history2->data;
		history2 = history2->down;
	}

	// a couple of variables for debug information
	double hypoth;
	// double second_best = 0;

	// time to iterate over all of the possible combinations. The avaliable
	// configurations are:
	// for all i in len_hist2 and all of j != k in len_hist1
	for (int i = 0; i < len_hist2; i++) {
		for (int j = 0; j < len_hist1; j++) {
			double i_j_dot = scatter_dir_dot(scatters2[i], scatters1[j]);
			vec3d* in = vec_sub(scatters2[i]->loc, scatters1[j]->loc);
			for (int k = 0; k < len_hist1; k++) {
				// can only do 3 point checks with j and k not being the same
				if (j != k) {
					// first check if the direction is physical
					if (((i_j_dot <= 0) && (i_j_dot != -1)) &&
							(scatter_dir_dot(scatters1[k], scatters1[j]) <= 0)) {
						vec3d* out = vec_sub(scatters1[j]->loc, scatters1[k]->loc);
						vec3d* gamma_cross = vec_cross(in, out);
						vec3d* gamma_cross_norm = vec_norm(gamma_cross);
						if (gamma_cross != NULL) {
							free(gamma_cross);
						}
						vec3d* ele_dir = vec_norm(scatters1[j]->dir);
						vec3d* plane = vec_cross(gamma_cross_norm, ele_dir);
						if (ele_dir != NULL)
							free(ele_dir);
						free(gamma_cross_norm);
						free(out);
						// the magnitude of plane is equal to the cos of the angle between the plane of the gamma
						// scattering and the direction of the electron

						if ((plane == NULL) || (vec_mag(plane) > 0.0)) {

							hypoth = expected_energy_b(scatters2[i], scatters1[j], scatters1[k], NULL);
							// check if the hypothesis is better than previous
							if (fabs(hypoth - ELECTRON_MASS) < fabs(best_find - ELECTRON_MASS)) {
								// the current hypthesis is better than the previous best
								if (GENERAL_DEBUG) {
									printf("multi_gamma_ele_iterator: new best scatter found:\n");
									printf("%f keV at ", hypoth);
									vec_print(scatters1[j]->loc, stdout);
									printf("\n");
								}
								// second_best = best_find;
								best_find = hypoth;
								best_scatter = scatters1[j];
								// predicted_vs_real[2 * (run_num - 1)] = j + 1;
								// predicted_vs_real[2 * (run_num - 1) + 1] = k + 1;
							}				
						}
						free(plane);
					}
				}
			}
			free(in);
		}
	}

	// if (first_scat_hypot == 0) {
	// 	first_scat_hypot = fabs(best_find - second_best);
	// } else {
	// 	second_scat_hypot = fabs(best_find - second_best);
	// }

	free(scatters1);
	free(scatters2);
	// done iterating, now have the best scatter found in the list
	if (fabs(best_find - ELECTRON_MASS) < (energy_percent * ELECTRON_MASS)) {
		// the result was within the energy cut
		if (GENERAL_DEBUG) {
			printf("multi_gamma_ele_iterator: best scatter solution found at %f kev at:\n", best_find);
			vec_print(best_scatter->loc, stdout);
			printf("\n");
		}
		return best_scatter;
	}
	// no result found within the energy cutoff
	return NULL;
}

/* build_array_no_i:
 * takes an array of scatters of length source_len and returns a copy of the
 * array without the element that was at location i in the source array.
 */
scatter** build_array_no_i(scatter** source, uint source_len, uint i) {
	if (source == NULL) {
		return NULL;
	}
	scatter** result = (scatter**)malloc(sizeof(scatter*) * (source_len - 1));
	for (int j = 0; j < source_len - 1; j++) {
		if (j < i) {
			result[j] = source[j];
		} else {
			result[j] = source[j+1];
		}
	}
	return result;
}

double recursive_search(double best, double current, double inc_eng, scatter* origin, scatter* loc, scatter** remaining, uint remain_count) {
	if ((loc == NULL) || (remaining == NULL) || (origin == NULL) ||(remain_count == 0)) {
		return current;
	}
	if (best < current) {
		// the best find is better than our current find, we are not doing well
		return current;
	}

	double OUT_OF_PLANE_WEIGHT = 2.; // how many sigma a 90 deg out of plane
	// scatter counts for 

	double energy_uncert;
	vec3d* in = vec_sub(origin->loc, loc->loc);

	double better_find = INFINITY;


	for (int i = 0; i < remain_count; i++) {
		// first calcuate the probability of our current deviation
		double expectation = expected_energy_b(origin, loc, remaining[i], &energy_uncert);
		if (energy_uncert < 0) {
			fprintf(stderr, "recursive_search: negative energy error found. Check given scatter errors\n");
		}
		double energy_error = (inc_eng - expectation) / energy_uncert; // error between
		double new_energy = inc_eng - loc->deposit; // could maybe be the value
		// of expectation - loc->deposit
		// the a->b gamma prediction and the value from earlier in the chain
		double step_error = add_quadrature(current, energy_error);

		// now find the error in electron direction plane
		if (loc->dir != NULL) {
			// only do this if there is a point. No electron direction then you
			// can't check if it is in the plane.

			vec3d* outgoing = vec_sub(loc->loc, remaining[i]->loc);
			// direction heading out from the current location to the next location
			vec3d* scatter_plane_normal = vec_cross(in, outgoing);
			// creates a normal to the plane of the scattering origin->loc->i
			vec3d* plane_norm_hat = vec_norm(scatter_plane_normal);
			// normalizes the plane perpendicular
			double out_of_plane = vec_dot(plane_norm_hat, loc->dir);
			// has the value of cosine of the angle between the plane normal and the
			// direction of the electron. This means that for perfect in the plane
			// behavior the value is 0, out of the plane +- 1.
			step_error = add_quadrature(step_error, abs(out_of_plane * OUT_OF_PLANE_WEIGHT));
			free(outgoing);
			free(scatter_plane_normal);
			free(plane_norm_hat);
		}
		// done calcuating what the current paths uncertanty is, lets combine it
		// with the current uncertanty to check if we are worse than the best
		// solution so far.
		double combined_error = add_quadrature(step_error, current);


		if (step_error < best) {
			// time to go one layer further
			scatter** new_array = build_array_no_i(remaining, remain_count, i);
			double below = recursive_search(best, combined_error, new_energy, loc, remaining[i], new_array, remain_count - 1);

			free(new_array);
			if (below < better_find) {
				better_find = below;
				if (below < best) {
					best = below;
				}
			}
		}
	}
	free(in);
	return better_find;
}

scatter* multi_gamma_stat_iteration(llist* history_near, llist* history_far, double sigma_per_scatter) {
	if ((history_near == NULL) || (history_far == NULL)) {
		return NULL;
	}
	// for each in history near, choose a history far and run the recursive search.
	// hang onto the best result after each point


	scatter* best_scatter = NULL;

	int len_hist_near = list_length(history_near);
	int len_hist_far = list_length(history_far);

	double best_find = sigma_per_scatter * len_hist_near;

	if (len_hist_near < 2) {
		// history1 is too short to run the recursion process.
		return NULL;
	}
	if (len_hist_far < 1) {
		// idk if it is possible to even get here, but we can't continue if we do
		return NULL;
	}
	llist* hist_near_bottom = list_tail(history_near);
	// time to turn the scattering list into an array for fast access
	scatter** scatters_near = (scatter**)malloc(len_hist_near * sizeof(scatter*));
	for (int i = 0; i < len_hist_near; i++)	{
		scatters_near[i] = (scatter*)hist_near_bottom->data;
		hist_near_bottom = hist_near_bottom->up;
	}

	llist* hist_far_bottom = list_tail(history_far);
	scatter** scatters_far = (scatter**)malloc(len_hist_far * sizeof(scatter*));
	for (int i = 0; i < len_hist_far; i++)	{
		scatters_far[i] = (scatter*)hist_far_bottom->data;
		hist_far_bottom = hist_far_bottom->up;
	}

	scatter_quicksort(scatters_near, 0, len_hist_near - 1);
	scatter_quicksort(scatters_far, 0, len_hist_far - 1);
	// sorts the scatters by deposit energy, now only keep the largest LARGEST
	scatter** scatters_near_short = (scatter**)malloc(LARGEST * sizeof(scatter*));
	scatter** scatters_far_short = (scatter**)malloc(LARGEST * sizeof(scatter*));
	for (int i = 0; i < LARGEST; i++) {
		if (i > len_hist_near) {
			scatters_near_short[i] = NULL;
		} else {
			scatters_near_short[i] = scatters_near[i];
		}
		if (i > len_hist_far) {
			scatters_far_short[i] = NULL;
		} else {
			scatters_far_short[i] = scatters_far[i];
		}
	}
	if (len_hist_near > LARGEST) {
		len_hist_near = LARGEST;
	}
	if (len_hist_far > LARGEST) {
		len_hist_far = LARGEST;
	}




	for (int j = 0; j < len_hist_near; j++) {

		scatter** new_array = build_array_no_i(scatters_near_short, len_hist_near, j);

		for (int i = 0; i < len_hist_far; i++) {

			double try = recursive_search(best_find, 0., 511., scatters_far_short[i], scatters_near_short[j], new_array, len_hist_near - 1);

			if (try < best_find) {
				best_scatter = scatters_near[j];
				best_find = try;
			}
		}
		free(new_array);
	}

	free(scatters_near);
	free(scatters_far);
	free(scatters_near_short);
	free(scatters_far_short);


	return best_scatter;

}

/*
 * closest_gamma:
 * finds the id of the closest gamma to the given location. The gamma is sourced
 * from a given history. The gamma is required to be from an annihilation
 * If the search fails returns zero
 */
int closest_gamma(llist* history, vec3d* target) {
	if ((history == NULL) || (target == NULL)) {
		return 0;
	}
	history = list_head(history);
	double distance = INFINITY;
	int best_id = 0;
	event* working_event = NULL;
	while (history != NULL) {
		// first, define our working event
		working_event = (event*)history->data;
		if ((working_event->particle == 22) && (working_event->orgin[0] == 'a')) {
			// it is a gamma from an annihilation
			double check_dist = vec_dist(target, working_event->location);
			if (distance > check_dist) {
				distance = check_dist;
				best_id = working_event->id;
			}
		}
		// move on to next event
		history = history->down;
	}
	return best_id;
}

/*
 * build_scatters
 * takes the list of interactions in the detector volume and makes a list of
 * energy deposition locations that are the scatter locations. This is the
 * current version of the "cluster finding" algorithm. It uses truth data to
 * find the location of the deposit. It also uses the real energy quantity
 * deposited at the location, not an estimate using the switchillator.
 * It returns the head of a new list of the scatter locations.
 */
llist* build_scatters(llist* detector_history, int id) {
	if (detector_history == NULL) {
		return NULL;
	}
	// looks to find the first instance of each electron showing up. Using the
	// location of that eele_dirith a scatter
	// of a gamma with the same identifier as id. If all is good the location is
	// added as a scatter with the electron's energy.

	// now also determines the direction of the electron. This is returned as
	// a simple vector pointing from the current place to the second location
	// of the electron. If the electron does not produce a second location then
	// the scatter is given a NULL value for the direction

	// put us at the start of the list
	detector_history = list_head(detector_history);
	// make the place for the list of scatters
	llist* scatter_list = NULL;
	// we will now iterate over the history looking for electrons. For each
	// electron we then find the associated gamma using closest_gamma. This is
	// a gamma match using only position. In theory the distance should be zero,
	// in pratice I don't trust Geant4 that much.
	llist* checked_electrons = NULL;

	while (detector_history != NULL) {
		int electron_checked = 0;
		// is it an electron, if not skip to next entry
		if (((event*)detector_history->data)->particle == 11) {
			// have we recorded this electron already?
			if (checked_electrons != NULL) {
				// list of checked electons is not empty
				llist* check = list_head(checked_electrons);
				while ((check != NULL) && (!electron_checked)) {
					if (*((int*)check->data) == ((event*)detector_history->data)->id) {
						// the current electron is in our list of handled variables
						electron_checked = 1;
					}
					check = check->down;
				}
			}
			// electron_checked now contains if we have previously looked at this
			// electron
			if (!electron_checked) {
				// if we haven't looked at this one yet, so is it compton based
				if (((event*)detector_history->data)->orgin[0] == 'c') {
					//  so which gamma is it from?
					int gamma_id = closest_gamma(detector_history, ((event*)detector_history->data)->location);
					if (gamma_id == id) {
						// time to add this electron to the scatter list
						vec3d* scatter_loc = vec_copy(((event*)detector_history->data)->location);
						event* cur_event = (event*)detector_history->data;
						event* next_event;
						if (detector_history->down != NULL) {
							next_event = (event*)detector_history->down->data;
						} else {
							next_event = NULL;
						}
						vec3d* ele_dir = NULL;
						if ((next_event != NULL) && (cur_event->id == next_event->id)) {
							ele_dir = vec_sub(next_event->location, ((event*)detector_history->data)->location);
							if ((fabs(ele_dir->x) < 0.0001) && (fabs(ele_dir->y) < 0.0001) && (fabs(ele_dir->z) < 0.0001)) {
								// no actual movement of the electron from which to find a path
								free(ele_dir);
								ele_dir = NULL;
							}
						} else {
							ele_dir = NULL;
						}
						scatter* add_scatter;
						if (ele_dir == NULL) {
							
							add_scatter = new_scatter(scatter_loc, NULL, cur_event->energy, cur_event->tof, sqrt(cur_event->energy), SPC_UNCERT, -1);
						} else {
							// normalize the electron direction
							vec3d* ele_dir_norm = vec_norm(ele_dir);
							free(ele_dir);
							add_scatter = new_scatter(scatter_loc, ele_dir_norm, cur_event->energy, cur_event->tof, sqrt(cur_event->energy), SPC_UNCERT, -1);
						}
						scatter_list = add_to_bottom(scatter_list, add_scatter);
					}
					// add the electron to the list of electrons we have checked
					int* electron_id = (int*)malloc(sizeof(int));
					*electron_id = ((event*)detector_history->data)->id;
					checked_electrons = add_to_bottom(checked_electrons, electron_id);
				}
			}
		}
		detector_history = detector_history->down;
	}

	fmap(checked_electrons, free_null);
	delete_list(checked_electrons);

	// go through and add the real n value for each scatter
	// the list will be in reverse order due to TOPAS/Geant4 call structure
	int length = list_length(scatter_list);
	llist* head = scatter_list;
	for (int i = 0; i < length; i++) {
		scatter_truth* true_info = ((scatter*)(head->data))->truth;
		if (true_info != NULL) {
			true_info->true_n = length - i;
		}
		head = head->down;
	}

	return scatter_list;
}

// we now have the infirsturcture for finding what scatter associated with a
// gamma is the first. Now we need to run this twice (once on each annihilation
// gamma) to get the two ends. Then from that we have the two ends and can
// run a modified version of first_scat_miss from the energy cut calculations

/*
 * find_endpoints
 * finds the ids of the first two gammas in the detector history. It is possible
 * for there to be more gammas than this due to intresting interactions, however
 * I am not looking to try and differentiate them currently. Once the first two
 * gamma ids are found, a call to build scatter is done on the first, then
 * scattering_iterator is run to find the expected first scatter. This is then
 * repeated for the other id. These two locations are then returned as a length
 * 2 array of the two scatter pointers. The line of responce is given by the
 * line between these two scatters
 */
scatter** find_endpoints(llist* detector_history, double energy_percent) {
	if (detector_history == NULL) {
		return NULL;
	}
	// find the id of the first gamma

	// first make sure we are at the top of the detector list
	detector_history = list_head(detector_history);
	llist* search_history = detector_history;
	int first_id = 0;
	int second_id = 0;

	while (search_history != NULL) {
		event* current_event = (event*)search_history->data;
		if ((current_event->particle == 22) && (current_event->id != first_id)) {
			// we have a gamma, and it isn't the same as the (possibly alreadly found)
			// first gamma
			if (!first_id) {
				first_id = current_event->id;
			} else {
				second_id = current_event->id;
				break;
			}
		}
		search_history = search_history->down;
	}
	if (second_id == 0) {
		// no second gamma was found, so no line of responce can be made.
		return NULL;
	}

	llist* scat_list = build_scatters(detector_history, first_id);
	if (scat_list == NULL) {
		return NULL;
	}
	scatter* endpoint = scattering_iterator(scat_list, energy_percent);
	// copy the scatter information over to a new structure so that the old list
	// can be freed
	if (endpoint == NULL) {
		return NULL;
	}
	scatter* first_endpoint = copy_scatter(endpoint);
	// free the old list of scatters:
	fmap(scat_list, delete_scatter);
	delete_list(scat_list);
	// repeat for the second scatter
	scat_list = build_scatters(detector_history, second_id);
	if (scat_list == NULL) {
		return NULL;
	}
	endpoint = scattering_iterator(scat_list, energy_percent);
	if (endpoint == NULL) {
		return NULL;
	}
	scatter* second_endpoint = copy_scatter(endpoint);
	fmap(scat_list, delete_scatter);
	delete_list(scat_list);

	// make an array of the two new scatters to be sent out of the function
	scatter** return_vals = (scatter**)malloc(2 * sizeof(scatter*));
	return_vals[0] = first_endpoint;
	return_vals[1] = second_endpoint;
	return return_vals;
}

/*
 * find_endpoints_2hist
 * finds the ids of the first two gammas in the detector history. This looks using
 * the list from the detector history, passing multi_gamma_iterator to look at
 * the two histories and find cases where hist2->hist1a->hist1b produce useable
 * energies.
 */
scatter** find_endpoints_2hist(llist* detector_history, double energy_percent) {
	if (detector_history == NULL) {
		return NULL;
	}
	// find the id of the first gamma

	// first make sure we are at the top of the detector list
	detector_history = list_head(detector_history);
	llist* search_history = detector_history;
	int first_id = 0;
	int second_id = 0;

	while (search_history != NULL) {
		event* current_event = (event*)search_history->data;
		if ((current_event->particle == 22) && (current_event->id != first_id)) {
			// we have a gamma, and it isn't the same as the (possibly alreadly found)
			// first gamma
			if (!first_id) {
				first_id = current_event->id;
			} else {
				second_id = current_event->id;
				break;
			}
		}
		search_history = search_history->down;
	}
	if (second_id == 0) {
		// no second gamma was found, so no line of responce can be made.
		return NULL;
	}

	llist* scat_list1 = build_scatters(detector_history, first_id);
	llist* scat_list2 = build_scatters(detector_history, second_id);
	if ((scat_list1 == NULL) || (scat_list2 == NULL)) {
		if (scat_list1 == NULL) {
			fmap(scat_list2, delete_scatter);
			delete_list(scat_list2);
		} else {
			fmap(scat_list1, delete_scatter);
			delete_list(scat_list1);
		}
		return NULL;
	}
	// run the actual finding of endpoint 1
	scatter* endpoint1 = multi_gamma_iterator(scat_list1, scat_list2, energy_percent);
	// now do the same with endpoint 2
	scatter* endpoint2 = multi_gamma_iterator(scat_list2, scat_list1, energy_percent);


	if ((endpoint1 == NULL) || (endpoint2 == NULL)) {
		fmap(scat_list1, delete_scatter);
		fmap(scat_list2, delete_scatter);
		delete_list(scat_list1);
		delete_list(scat_list2);
		return NULL;
	}

	if (GENERAL_DEBUG) {
		printf("find_endpoints_2hist: scatters found at:\n");
		vec_print(endpoint1->loc, stdout);
		printf("\n");
		vec_print(endpoint2->loc, stdout);
		printf("\n");
	}

	// copy the two endpoints to new scatter structures (allows freeing of
	// scatter lists)

	scatter *first_endpoint = copy_scatter(endpoint1);
	scatter *second_endpoint = copy_scatter(endpoint2);

	// free the old list of scatters:
	fmap(scat_list1, delete_scatter);
	fmap(scat_list2, delete_scatter);
	delete_list(scat_list1);
	delete_list(scat_list2);
	

	// make an array of the two new scatters to be sent out of the function
	scatter** return_vals = (scatter**)malloc(2 * sizeof(scatter*));
	return_vals[0] = first_endpoint;
	return_vals[1] = second_endpoint;
	if (GENERAL_DEBUG) {
		printf("find_endpoints_2hist: returning scatters found at:\n");
		vec_print(return_vals[0]->loc, stdout);
		printf("\n");
		vec_print(return_vals[1]->loc, stdout);
		// vec_print(first_endpoint->loc, stdout);
		// printf("\n");
		// vec_print(second_endpoint->loc, stdout);
		
		printf("\n");
	}
	return return_vals;
}

/*
 * find_endpoints_ele_dir
 * finds the ids of the first two gammas in the detector history. This looks using
 * the list from the detector history, passing multi_gamma_iterator to look at
 * the two histories and find cases where hist2->hist1a->hist1b produce useable
 * energies. Also takes into account the electron direction to reject
 * non-physical solutions where the electron is in the same direction as the
 * vector from the scatter to the origin or resulting location
 */
scatter** find_endpoints_ele_dir(llist* detector_history, double energy_percent) {
	if (detector_history == NULL) {
		return NULL;
	}
	// find the id of the first gamma

	// first make sure we are at the top of the detector list
	detector_history = list_head(detector_history);
	llist* search_history = detector_history;
	int first_id = 0;
	int second_id = 0;

	while (search_history != NULL) {
		event* current_event = (event*)search_history->data;
		if ((current_event->particle == 22) && (current_event->id != first_id)) {
			// we have a gamma, and it isn't the same as the (possibly alreadly found)
			// first gamma
			if (!first_id) {
				first_id = current_event->id;
			} else {
				second_id = current_event->id;
				break;
			}
		}
		search_history = search_history->down;
	}
	if (second_id == 0) {
		// no second gamma was found, so no line of responce can be made.
		return NULL;
	}

	llist* scat_list1 = build_scatters(detector_history, first_id);
	llist* scat_list2 = build_scatters(detector_history, second_id);
	if ((scat_list1 == NULL) || (scat_list2 == NULL)) {
		if (scat_list1 == NULL) {
			fmap(scat_list2, delete_scatter);
			delete_list(scat_list2);
		} else {
			fmap(scat_list1, delete_scatter);
			delete_list(scat_list1);
		}
		return NULL;
	}
	// run the actual finding of endpoint 1
	scatter* endpoint1 = multi_gamma_ele_iterator(scat_list1, scat_list2, energy_percent);
	// now do the same with endpoint 2
	scatter* endpoint2 = multi_gamma_ele_iterator(scat_list2, scat_list1, energy_percent);


	if ((endpoint1 == NULL) || (endpoint2 == NULL)) {
		fmap(scat_list1, delete_scatter);
		fmap(scat_list2, delete_scatter);
		delete_list(scat_list1);
		delete_list(scat_list2);
		return NULL;
	}

	if (GENERAL_DEBUG) {
		printf("find_endpoints_ele_dir: scatters found at:\n");
		vec_print(endpoint1->loc, stdout);
		printf("\n");
		vec_print(endpoint2->loc, stdout);
		printf("\n");
	}

	// copy the two endpoints to new scatter structures (allows freeing of
	// scatter lists)

	scatter *first_endpoint = copy_scatter(endpoint1);
	scatter *second_endpoint = copy_scatter(endpoint2);

	// free the old list of scatters:
	fmap(scat_list1, delete_scatter);
	fmap(scat_list2, delete_scatter);
	delete_list(scat_list1);
	delete_list(scat_list2);
	

	// make an array of the two new scatters to be sent out of the function
	scatter** return_vals = (scatter**)malloc(2 * sizeof(scatter*));
	return_vals[0] = first_endpoint;
	return_vals[1] = second_endpoint;
	if (GENERAL_DEBUG) {
		printf("find_endpoints_ele_dir: returning scatters found at:\n");
		vec_print(return_vals[0]->loc, stdout);
		printf("\n");
		vec_print(return_vals[1]->loc, stdout);
		// vec_print(first_endpoint->loc, stdout);
		// printf("\n");
		// vec_print(second_endpoint->loc, stdout);
		
		printf("\n");
	}
	return return_vals;
}

/* 
 * find_endpoints_stat
 * finds the predicted endpoints of the first two gammas in the detector history
 * using a statistical determination method. For the first step a scatter is
 * chosen from each gamma. Then recursivly a single scatter is chosen from the
 * list of remaining scatters being seached, and the probability of the deviation
 * from physical that is found is determined. If this does not exceed the current
 * best find the process repeats until the list is empty. The recursive portion
 * is completed using the recursive_search function. The energy_percent value
 * provides a minimum probability that must be cleared for the result to be
 * passed out.
 */
scatter** find_endpoints_stat(llist* detector_history, double sigma_per_scatter) {
	if (detector_history == NULL) {
		return NULL;
	}
	// find the id of the first gamma

	// first make sure we are at the top of the detector list
	detector_history = list_head(detector_history);
	llist* search_history = detector_history;
	int first_id = 0;
	int second_id = 0;

	while (search_history != NULL) {
		event* current_event = (event*)search_history->data;
		if ((current_event->particle == 22) && (current_event->id != first_id)) {
			// we have a gamma, and it isn't the same as the (possibly alreadly found)
			// first gamma
			if (!first_id) {
				first_id = current_event->id;
			} else {
				second_id = current_event->id;
				break;
			}
		}
		search_history = search_history->down;
	}
	if (second_id == 0) {
		// no second gamma was found, so no line of responce can be made.
		return NULL;
	}

	llist* scat_list1 = build_scatters(detector_history, first_id);
	llist* scat_list2 = build_scatters(detector_history, second_id);
	if ((scat_list1 == NULL) || (scat_list2 == NULL)) {
		if (scat_list1 == NULL) {
			fmap(scat_list2, delete_scatter);
			delete_list(scat_list2);
		} else {
			fmap(scat_list1, delete_scatter);
			delete_list(scat_list1);
		}
		return NULL;
	}
	// run the actual finding of endpoint 1
	scatter* endpoint1 = multi_gamma_stat_iteration(scat_list1, scat_list2, sigma_per_scatter);
	// now do the same with endpoint 2
	scatter* endpoint2 = multi_gamma_stat_iteration(scat_list2, scat_list1, sigma_per_scatter);


	if ((endpoint1 == NULL) || (endpoint2 == NULL)) {
		fmap(scat_list1, delete_scatter);
		fmap(scat_list2, delete_scatter);
		delete_list(scat_list1);
		delete_list(scat_list2);
		return NULL;
	}

	if (GENERAL_DEBUG) {
		printf("find_endpoints_ele_dir: scatters found at:\n");
		vec_print(endpoint1->loc, stdout);
		printf("\n");
		vec_print(endpoint2->loc, stdout);
		printf("\n");
	}

	// copy the two endpoints to new scatter structures (allows freeing of
	// scatter lists)

	scatter *first_endpoint = copy_scatter(endpoint1);
	scatter *second_endpoint = copy_scatter(endpoint2);

	// free the old list of scatters:
	fmap(scat_list1, delete_scatter);
	fmap(scat_list2, delete_scatter);
	delete_list(scat_list1);
	delete_list(scat_list2);
	

	// make an array of the two new scatters to be sent out of the function
	scatter** return_vals = (scatter**)malloc(2 * sizeof(scatter*));
	return_vals[0] = first_endpoint;
	return_vals[1] = second_endpoint;
	if (GENERAL_DEBUG) {
		printf("find_endpoints_ele_dir: returning scatters found at:\n");
		vec_print(return_vals[0]->loc, stdout);
		printf("\n");
		vec_print(return_vals[1]->loc, stdout);
		// vec_print(first_endpoint->loc, stdout);
		// printf("\n");
		// vec_print(second_endpoint->loc, stdout);
		
		printf("\n");
	}
	return return_vals;
}


/* 
 * first_scat_miss:
 * Takes the two endpoints and a 3-vector defining the location of the
 * annihilation and finds the distance between the line-of-responce and the
 * annihilation. If there is a problem (such as not having two endpoints) it
 * returns -1 as a reject value.
 */
double first_scat_miss(scatter** endpoints, vec3d* annh_loc) {
	if ((endpoints == NULL) || (annh_loc == NULL)) {
		return -1.;
	}
	if ((endpoints[0] == NULL) || (endpoints[1] == NULL)) {
		return -1.;
	}
	// find the scattering point of the first scatter

	vec3d* first_spot = vec_copy(endpoints[0]->loc);

	vec3d* second_spot = vec_copy(endpoints[1]->loc);

	return line_to_dot_dist(first_spot, second_spot, annh_loc);
}

/* 
 * looks to see if the given event is from the annihilation.
 * If it is, it returns the particle id, otherwise it returns a zero
 */
int find_annih_gamma(event* item) {
	if (item == NULL) {
		fprintf(stderr, "find_annih_gamma passed NULL pointer\n");
		return 0;
	}
	if (item->orgin[0] == 'a') {
		return item->id;
	}
	return 0;
}

/*
 * create_lor
 * Takes two events and returns a LOR based on the two locations. The LOR is
 * defined by being between the two endpoints.
 * Method:
 * 		To find the centerpoint it finds the physical center between the two
 * 		locations, the displaces that along the line based on the difference
 * 		between the times of the two endpoints. The displacement is by 
 * 		c * (delta T) in the direction of the lower time.
 * 		
 * 		The length uncertanty is based on the time uncertanty of the two
 * 		scatters and the space uncertanty
 * 		The width uncertanty is based on the space uncertanty of the two
 * 		scatters
 */
lor* create_lor(scatter* a, scatter* b) {
	if ((a == NULL) || (b == NULL)) {
		return NULL;
	}
	vec3d* center_subtraction = vec_sub(a->loc, b->loc);
	vec3d* center_half = vec_scaler(center_subtraction, 0.5);
	vec3d* geometric_center = vec_add(b->loc, center_half);
	vec3d* ba_unit = vec_norm(center_half);
	double time_delta = b->time - a->time;
	vec3d* displacement = vec_scaler(ba_unit, SPD_LGHT * time_delta * 0.5);
	lor* new = (lor*)malloc(sizeof(lor));
	new->center = vec_add(geometric_center, displacement);
	new->dir = ba_unit;
	double a_space = a->space_uncert;
	double b_space = b->space_uncert;
	double a_time  = a->time_uncert;
	double b_time  = b->time_uncert;
	if (a_space < 0) {
		a_space = .1;
	}
	if (b_space < 0) {
		b_space = 0.1;
	}
	if (a_time < 0) {
		a_time = 5.;
	}
	if (b_time < 0) {
		b_time = 5.;
	}

	new->transverse_uncert = sqrt(a_space * a_space + b_space * b_space);
	new->long_uncert = sqrt(a_space * a_space + b_space * b_space + 
							a_time * a_time + b_time * b_time);
	free(center_subtraction);
	free(center_half);
	free(geometric_center);
	free(displacement);
	return new;
}

void print_lor(FILE* output, lor* lor) {
	fprintf(output, "%f, %f, %f,", lor->center->x, lor->center->y, lor->center->z);
	fprintf(output,  " %f, %f, %f,", lor->dir->x, lor->dir->y, lor->dir->z);
	fprintf(output, " %f, %f", lor->long_uncert, lor->transverse_uncert);
}

struct _source_union {
	llist* history;
	double cutoff;
	FILE* output;
};

// a wrapper of the statistical inverse kinematics algoritem that allows for
// multithreading. Locks the sending to fprintf.
void* wrapper_inv_kin_stat(void* a) {
	struct _source_union *source = (struct _source_union *)a;
	llist* history = source->history;
	FILE* lor_output = source->output;
	scatter** endpoints = find_endpoints_stat(history, source->cutoff);

	if (endpoints == NULL) {

	} else {
		lor* result = create_lor(endpoints[0], endpoints[1]);
		pthread_mutex_lock(&file_lock);
		fprintf(lor_output, "%i, ", ((event*)(history->data))->number);
		print_lor(lor_output, result);
		fprintf(lor_output, "\n");
		pthread_mutex_unlock(&file_lock);

		free_lor(result);

		delete_scatter(endpoints[0]);
		delete_scatter(endpoints[1]);
		free(endpoints);

	}
	fmap(history, delete_event);
	delete_list(history);
	free(source);
	return 0;
}


int main(int argc, char **argv) {
	// quick reminder that argc is the number of arguments,
	// argv points to the strings of the arguments

	// we are looking for an input histories file, an output file name,
	// and an inner radius of the detector this may later be changed to
	// input histories, input data file, output file name
	if (argc != 4) {
		printf("Unable to run. Expected 3 arguments got %i.\n", argc - 1);
		printf("Expects an input detector volume ");
		printf("history file, output file name, and");
		printf(" energy cutoff percent.\n");
		return 1;
	}
	// FILE* in_histories = fopen(argv[1], "r");
	// if (in_histories == NULL) {
	// 	printf("Unable to open history file\n");
	// 	return 1;
	// }
	FILE* in_det_histories = fopen(argv[1], "r");
	if (in_det_histories == NULL) {
		printf("Unable to open detector history file\n");
		return 1;
	}
	char* lor_file = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
	strcpy(lor_file, argv[2]);
	lor_file = strcat(lor_file, ".lor");
	FILE* lor_output = fopen(lor_file, "w");
	if (lor_output == NULL) {
		printf("Unable to open output file for writing\n");
		return 1;
	}

	energy_cutoff = strtod(argv[3], NULL);

	if (!test_expected_energy()) {
		fprintf(stderr, "tests failed, exiting\n");
		return 1;
	}


	llist *in_det_hist = load_historyb(in_det_histories, read_line);
	uint run_num = 0;
	first_scat_hypot = 0;
	second_scat_hypot = 0;

	if (pthread_mutex_init(&file_lock, NULL)) {
		fprintf(stderr, "pthread: unable to make volume lock, exiting.\n");
		return(1);
	}

	for (int i = 0; i < MAX_THREAD_CALLS; i++) {
		// set the array of threads to be empty
		tid[i] = -1;
	}
	int cur_thread = 0;


	// begin the primary loop over all histories
	while (in_det_hist != NULL) {
		// print an update to how far we have made it
		if ((run_num / 10000) * 10000 == run_num) {
			printf("run number: %u\n", run_num);
		}
		// if (tid[cur_thread] != -1) {
		// 	pthread_join(tid[cur_thread], NULL);
		// 	tid[cur_thread] = -1;
		// }

		// lets make this multithreaded
		// struct _source_union *arguments = (struct _source_union *)malloc(sizeof(struct _source_union));
		// arguments->history = in_det_hist;
		// arguments->cutoff = energy_cutoff;
		// arguments->output = lor_output;
		// pthread_create(&tid[cur_thread], NULL, wrapper_inv_kin_stat, arguments);

		scatter** endpoints = find_endpoints_stat(in_det_hist, energy_cutoff);

		if (endpoints == NULL) {

		} else {
			// create the LOR
			lor* result = create_lor(endpoints[0], endpoints[1]);
			fprintf(lor_output, "%i, ", ((event*)(in_det_hist->data))->number);
			print_lor(lor_output, result);
			fprintf(lor_output, "\n");
		}
		if (endpoints != NULL) {
			delete_scatter(endpoints[0]);
			delete_scatter(endpoints[1]);
			free(endpoints);
		}
		fmap(in_det_hist, delete_event);
		delete_list(in_det_hist);
		in_det_hist = load_historyb(in_det_histories, read_line);
		run_num++;

		// if (cur_thread + 1 < MAX_THREAD_CALLS) {
		// 	cur_thread++;
		// } else {
		// 	cur_thread = 0;
		// }
		// scatter** endpoints = find_endpoints_ele_dir(in_det_hist, energy_cutoff);
		// scatter** endpoints = find_endpoints_stat(in_det_hist, energy_cutoff);

		// if (endpoints == NULL) {

		// } else {
			
		// 	lor* result = create_lor(endpoints[0], endpoints[1]);
		// 	fprintf(lor_output, "%i, ", ((event*)(in_det_hist->data))->number);
		// 	print_lor(lor_output, result);
		// 	fprintf(lor_output, "\n");
		// 	free_lor(result);

		// 	delete_scatter(endpoints[0]);
		// 	delete_scatter(endpoints[1]);
		// 	free(endpoints);
		// }
		// fmap(in_det_hist, delete_event);
		// delete_list(in_det_hist);
		// in_det_hist = load_historyb(in_det_histories, read_line);
		// run_num++;

	}


	fclose(in_det_histories);
	fclose(lor_output);
	return 0;
}

