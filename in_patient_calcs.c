#include "truth_assign.h"
#include <stdio.h>
#include <string.h>
#include "llist.h"
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"

#define ENG_RNG 0.001
#define COMP_INT 1667457891
#define ELECTRON_MASS 510.999


  
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

// creates a new scatter structure and fills it
scatter* new_scatter_old(vec3* vector, double deposited) {
	scatter* new = (scatter*)malloc(sizeof(scatter));
	new->deposit = deposited;
	new->loc = vector;
	new->eng_uncert = -1.;
	new->space_uncert = -1.;
	return new;
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

// frees scatter malloc
void* delete_scatter(void* in) {
	if (in == NULL) {
		return NULL;
	}
	free(((scatter*)in)->loc);
	free(in);
	return NULL;
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
	// make a new copy of the event in previous event for storing
	// means it will continue pointing right as otherwise it can
	// point to the history that got freed
	previous_event = duplicate_event(previous_event);
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
vec3* find_annihilation_point(llist *history) {
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
double line_to_dot_dist(vec3* start, vec3* end, vec3* point) {
	vec3* num_first_term = vec_sub(start, end);
	vec3* num_sec_term = vec_sub(start, point);
	// vec3* denom_vec = vec_sub(end, start);
	double numerator = vec_mag(vec_cross(num_first_term, num_sec_term));
	double denomenator = vec_mag(num_first_term);
	free(num_first_term);
	free(num_sec_term);
	// free(denom_vec);
	return numerator / denomenator;
}


/* 
 * expected_energy_b
 * takes 3 scatters and solves the kinematics assuming that the gamma goes from
 * a->b->c, solving for the energy of the gamma from a->b. To do this it
 * calculates the angle <ABC, then uses the deposited energy at B and this angle
 * to determine the energy of the incoming gamma. 
 */
 double expected_energy_b(scatter* a, scatter* b, scatter* c) {
	 if ((a == NULL) || (b == NULL) || (c == NULL)) {
		 return 1.;
	 }

	// first calculate the angle at b

	// get the vector from b->a
	vec3* ab = vec_sub(b->loc, a->loc);
	// get the vector from b->c
	vec3* bc = vec_sub(c->loc, b->loc);
	// calculate the angle itself
	double theta = vec_angle(ab, bc);

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
	double gamma_to_b_e = expected_energy_b(a, b, c);

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

	vec3* point_a = three_vec(0., 0., 0.);
	vec3* point_b = three_vec(0., 3., 0);
	vec3* point_c = three_vec(0., 2., 1.73205);
	double deposit_a = 127.405;
	double deposit_b = 203.1654;
	scatter* scatter_a = new_scatter_old(point_a, deposit_a);
	scatter* scatter_b = new_scatter_old(point_b, deposit_b);
	scatter* scatter_c = new_scatter_old(point_c, 10.);
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
scatter* scattering_iterator(llist* history, double energy_percent) {
	if (history == NULL) {
		fprintf(stderr, "scattering_iterator: no history given\n");
		return NULL;
	}

	// the current best find energy
	double best_find = 0;
	// the current best find scatter
	scatter* best_scatter = NULL;
	// length of the scatter list
	int list_len = list_length(history);

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
		locations[i] = (scatter*)history->data;
		history = history->down;
	}
	


	double hypoth;
	double second_best;
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
						second_best = best_find;
						best_find = hypoth;
						best_scatter = locations[first];
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

	scatter* best_scatter;

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
	double second_best;

	// time to iterate over all of the possible combinations. The avaliable
	// configurations are:
	// for all i in len_hist2 and all of j != k in len_hist1
	for (int i = 0; i < len_hist2; i++) {
		for (int j = 0; j < len_hist1; j++) {
			for (int k = 0; k < len_hist1; k++) {
				// can only do 3 point checks with j and k not being the same
				if (j != k) {
					hypoth = expected_energy_b(scatters2[i], scatters1[j], scatters1[k]);
					// check if the hypothesis is better than previous
					if (fabs(hypoth - ELECTRON_MASS) < fabs(best_find - ELECTRON_MASS)) {
						// the current hypthesis is better than the previous best
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
		return best_scatter;
	}
	// no result found within the energy cutoff
	return NULL;
}

/*
 * closest_gamma:
 * finds the id of the closest gamma to the given location. The gamma is sourced
 * from a given history. The gamma is required to be from an annihilation
 * If the search fails returns zero
 */
int closest_gamma(llist* history, vec3* target) {
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
	// location of that electron it checks to see if it lines up with a scatter
	// of a gamma with the same identifier as id. If all is good the location is
	// added as a scatter with the electron's energy.

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
						vec3* scatter_loc = vec_copy(((event*)detector_history->data)->location);
						scatter_list = add_to_bottom(scatter_list, new_scatter_old(scatter_loc, ((event*)detector_history->data)->energy));
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
	scatter* first_endpoint = new_scatter_old(vec_copy(endpoint->loc), endpoint->deposit);
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
	scatter* second_endpoint = new_scatter_old(vec_copy(endpoint->loc), endpoint->deposit);
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

	// copy the two endpoints to new scatter structures (allows freeing of
	// scatter lists)

	scatter *first_endpoint = new_scatter_old(endpoint1->loc, endpoint1->deposit);
	scatter *second_endpoint = new_scatter_old(endpoint2->loc, endpoint2->deposit);

	// free the old list of scatters:
	fmap(scat_list1, delete_scatter);
	fmap(scat_list2, delete_scatter);
	delete_list(scat_list1);
	delete_list(scat_list2);
	

	// make an array of the two new scatters to be sent out of the function
	scatter** return_vals = (scatter**)malloc(2 * sizeof(scatter*));
	return_vals[0] = first_endpoint;
	return_vals[1] = second_endpoint;
	return return_vals;
}

/* 
 * first_scat_miss:
 * Takes the two endpoints and a 3-vector defining the location of the
 * annihilation and finds the distance between the line-of-responce and the
 * annihilation. If there is a problem (such as not having two endpoints) it
 * returns -1 as a reject value.
 */
double first_scat_miss(scatter** endpoints, vec3* annh_loc) {
	if ((endpoints == NULL) || (annh_loc == NULL)) {
		return -1.;
	}
	if ((endpoints[0] == NULL) || (endpoints[1] == NULL)) {
		return -1.;
	}
	// find the scattering point of the first scatter

	vec3* first_spot = vec_copy(endpoints[0]->loc);

	vec3* second_spot = vec_copy(endpoints[1]->loc);

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


int main(int argc, char **argv) {
	// quick reminder that argc is the number of arguments,
	// argv points to the strings of the arguments

	// we are looking for an input histories file, an output file name,
	// and an inner radius of the detector this may later be changed to
	// input histories, input data file, output file name
	if (argc != 7) {
		printf("Unable to run. Expected 6 arguments got %i.\n", argc - 1);
		printf("Expects an input full history file, input detector volume ");
		printf("history file, output file name, detector");
		printf(" inner radius, detector half height and energy cutoff percent.\nDetector sizes");
		printf(" are in the same units as the input history file's distances.\n");
		return 1;
	}
	FILE* in_histories = fopen(argv[1], "r");
	if (in_histories == NULL) {
		printf("Unable to open history file\n");
		return 1;
	}
	FILE* in_det_histories = fopen(argv[2], "r");
	if (in_det_histories == NULL) {
		printf("Unable to open detector history file\n");
		return 1;
	}
	FILE* out_in_patient = fopen(argv[3], "w");
	if (out_in_patient == NULL) {
		printf("Unable to open output file for writing\n");
		return 1;
	}
	double in_patient_distance = strtod(argv[4], NULL);
	double detector_height = strtod(argv[5], NULL);
	energy_cutoff = strtod(argv[6], NULL);
	if ((in_patient_distance <= 0) || (detector_height <= 0) || (energy_cutoff <= 0)) {
		printf("Size and energy dimensions must be greater than zero\n");
		return 1;
	}
	if (!test_expected_energy()) {
		fprintf(stderr, "tests failed, exiting\n");
		return 1;
	}

	fprintf(out_in_patient, "history number, in patient scatter occurance, ");
	fprintf(out_in_patient, "first algorithem miss distance, second algo miss dist, ");
	fprintf(out_in_patient, "delta first hypot (2ed algo), delta sec hypot (2ed algo)\n");

	// current loop version, just for testing
	llist *history = load_history(in_histories, read_line);
	llist *in_det_hist = load_historyb(in_det_histories, read_line);
	// some extra vars for calculating percentage of in paitent scatters
	uint table[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	uint table_total = 0;
	// uint numerator = 0;
	// uint denominator = 0;
	uint hist_num = 0;
	first_scat_hypot = 0;
	second_scat_hypot = 0;
	while (history != NULL) {
		// print an update to how far we have made it
		if ((hist_num / 10000) * 10000 == hist_num) {
			printf("history number: %u\n", hist_num);
		}

		// find the first annihilation gamma
		int part_id = 0; // particle id
		int in_patient_occurance = 0;
		first_scat_hypot = 0;
		second_scat_hypot = 0;
		llist* curr_loc = list_head(history);
		while (!part_id) {
			part_id = find_annih_gamma((event*)curr_loc->data);
			curr_loc = curr_loc->down;
		}
		// does the first annhilation gamma scatter, and how much
		int scatters1 = in_patient(history, in_patient_distance, detector_height, part_id);
		in_patient_occurance = in_patient_occurance || scatters1;
		// find the second gamma 
		if (scatters1 < 3) {
			int previous = part_id;
			while ((part_id == previous) || (!part_id)) {
				part_id = find_annih_gamma((event*)curr_loc->data);
				curr_loc = curr_loc->down;
			}
			// does the second annhilation gamma scatter, and how much
			int scatters2 = in_patient(history, in_patient_distance, detector_height, part_id);
			in_patient_occurance = in_patient_occurance || scatters2;

			if ((scatters1 < 3) && (scatters2 < 3)) {
				// add one to the table of fates if not beyond measurements
				table[scatters1][scatters2]++;
				table_total++;
			}
		}
		// fmap(history, print_event);
		if (in_patient_occurance) {
			fprintf(out_in_patient, "%i, 1, ", ((event*)(history->data))->number);
		} else {
			fprintf(out_in_patient, "%i, 0, ", ((event*)history->data)->number);
		}
		// done determining in patient scattering

		// only can run the endpoint finding if main history and detector history
		// from the same history
		if (((event*)in_det_hist->data)->number == ((event*)history->data)->number) {
			// first we need to find the location of the endpoint scatters

			// first run on old code
			scatter** endpoints = find_endpoints(in_det_hist, energy_cutoff);

			if (endpoints == NULL) {
				fprintf(out_in_patient, "%f, ", -1.0);
			} else {
					// now we need to find the distance by which the endpoints miss
				vec3* annh_loc = find_annihilation_point(history);
				double miss_dist = first_scat_miss(endpoints, annh_loc);
				if (endpoints != NULL) {
					delete_scatter(endpoints[0]);
					delete_scatter(endpoints[1]);
					free(endpoints);
				}

				fprintf(out_in_patient, "%f, ", miss_dist);
			}
			first_scat_hypot = 0;
			second_scat_hypot = 0;

			endpoints = find_endpoints_2hist(in_det_hist, energy_cutoff);

			if (endpoints == NULL) {
				fprintf(out_in_patient, "%f, %f, %f\n", -1.0, -1.0, -1.0);
			} else {
					// now we need to find the distance by which the endpoints miss
				vec3* annh_loc = find_annihilation_point(history);
				double miss_dist = first_scat_miss(endpoints, annh_loc);
				if (endpoints != NULL) {
					// delete_scatter(endpoints[0]);
					// delete_scatter(endpoints[1]);
					free(endpoints);
				}

				fprintf(out_in_patient, "%f, ", miss_dist);
				// fprintf(out_in_patient, "%f, %f\n", first_scat_hypot, second_scat_hypot);
				vec_print(endpoints[0]->loc, out_in_patient);
				vec_print(endpoints[1]->loc, out_in_patient);
				fprintf(out_in_patient, "\n");
			}

		} else {
			fprintf(out_in_patient, "%f, %f, %f, %f\n", -1.0, -1.0, -1.0, -1.0);
		}


		fmap(history, delete_event);
		delete_list(history);
		history = load_history(in_histories, read_line);
		hist_num++;
		if (history != NULL) {		
			if (((event*)history->data)->number > ((event*)in_det_hist->data)->number) {
				fmap(in_det_hist, delete_event);
				delete_list(in_det_hist);
				in_det_hist = load_historyb(in_det_histories, read_line);
			}
		}

	}
	// if (denominator != 0) {
	// 	printf("Percentage of in paitent scatters: ");
	// 	float percent = (float)numerator / (float)denominator;
	// 	printf("%f\n", percent);
	// }

	printf("Total events: %u\nFates:\n", hist_num);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			printf("%4.4f		", (double)table[i][j]/(double)hist_num);
		}
		printf("\n");
	}


	fclose(in_histories);
	fclose(out_in_patient);
	return 0;
}

