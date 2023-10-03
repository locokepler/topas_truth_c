#include "truth_assign.h"
#include <stdio.h>
#include <string.h>
#include "llist.h"
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"

#define ENG_RNG 0.001
#define COMP_INT 1667457891
#define ELECTRON_MASS 510.999
#define SPD_LGHT 29.98
#define PI 3.141592653589
#define MIN_SCAT_ENG 10.0 // minumum energy for a scatter to be observed
#define E_TRIGGER 20.0 // energy in keV to hit trigger 

#define TIME_UNCERT_CM 3.82
#define SPC_UNCERT_PLANE 0.3 // uncertainty in the plane of the bore
#define SPC_UNCERT_RAD 1.0  // uncertainty radial to the bore
#define UNCERT_REP 12 // repetitions of addition to create random values
#define P_PER_KEV 30 // photons per keV. Efficiencys have already been applied
// about 50 for NaI(Tl)
// 30 for LYSO
#define BORE_RADIUS 45 // measured in cm
#define SCANNER_THICKNESS 9 // measured in cm
#define CRYSTAL_IN_PHI 1024 // these values roughly match uEXPLORER
#define CRYSTAL_HEIGHT 1.8
#define CRYSTAL_Z_THICKNESS 0.276


#define READ_DEBUG 0
#define GENERAL_DEBUG 0
#define TREE_DEBUG 0
#define GRAPHVIZ_DEBUG 0
#define CLUSTER_DEBUG 0
FILE* debug_graphs = NULL;
uint graph_id;


uint hist_num;
uint true_counts = 0;
uint scatter_counts = 0;

int correct_clustering = 0;
int total_clustering = 0;

typedef struct _crystal {
	double energy;
	double time;
	int index;
	uint true_n;
	int true_gamma;
	vec3d pos;
} crystal;
  
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

// reads an event from a binary file
event* read_line_binary(FILE* source) {

	// begin to fill the event struct with information from .phsp
	uint numb;
	double energy;
	double deposit;
	float x;
	float y;
	float z;
	float tof;
	int particle;
	int count;

	int worked = 0;

	worked += fread(&numb, sizeof(uint), 1, source);
	worked += fread(&energy, sizeof(double), 1, source);
	worked += fread(&deposit, sizeof(double), 1, source);
	worked += fread(&x, sizeof(float), 1, source);
	worked += fread(&y, sizeof(float), 1, source);
	worked += fread(&z, sizeof(float), 1, source);
	worked += fread(&tof, sizeof(float), 1, source);
	worked += fread(&particle, sizeof(int), 1, source);
	worked += fread(&count, sizeof(int), 1, source);

	if (worked != 9) {
		// something went wrong with the read, don't pass bad information
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
	new_event->location 	= three_vec((double)x,(double)y,(double)z);
	new_event->tof 			= (double)tof;
	new_event->particle 	= particle;
	new_event->id			= count;
	new_event->orgin[0]		= (char)0;

	// if (READ_DEBUG) {
	// 	print_event((void*)new_event);
	// }

	return new_event;
}

// creates a new scatter structure and fills it
scatter* new_scatter_old(vec3d vector, double deposited, double time) {
	scatter* new = (scatter*)malloc(sizeof(scatter));
	new->deposit = deposited;
	new->loc = vector;
	new->dir = three_vec(0,0,0);
	new->has_dir = 0;
	new->time = time;
	new->eng_uncert = -1;
	new->space_uncert = -1;
	new->time_uncert = -1;
	new->truth = NULL;
	return new;
}

scatter* new_scatter(vec3d vector, vec3d dir, char has_dir, double deposit, double time, double eng_uncert, double space_uncert, double time_uncert) {
	scatter* new = (scatter*)malloc(sizeof(scatter));
	new->deposit = deposit;
	new->eng_uncert = eng_uncert;
	new->space_uncert = space_uncert;
	new->loc = vector;
	new->dir = dir;
	new->has_dir = has_dir;
	new->time = time;
	new->time_uncert = time_uncert;
	new->truth = NULL;
	return new;
}

scatter* copy_scatter(scatter* a) {
	scatter* new = new_scatter(a->loc, a->dir, a->has_dir, a->deposit, a->time, a->eng_uncert, a->space_uncert, a->time_uncert);
	if (a->truth != NULL) {
		new->truth = (scatter_truth*)malloc(sizeof(scatter_truth));
		new->truth->true_eng  = a->truth->true_eng;
		new->truth->true_n    = a->truth->true_n;
		new->truth->true_time = a->truth->true_time;
	}
	return new;
}

void print_scatter(scatter* a) {
	if (a == NULL) {
		return;
	}
	printf("scatter: E = %lf +- %lf, pos: ", a->deposit, a->eng_uncert);
	vec_print(a->loc, stdout);
	printf(", t = %lf +- %lf, gamma %i", a->time, a->time_uncert, (int)(a->has_dir));
	if (a->truth != NULL) {
		printf(", N = %X", a->truth->true_n);
	}
	printf("\n");
}

void* mappable_print_scatter(void* a) {
	print_scatter((scatter*)a);
	return a;
}

void* delete_scatter(void* in) {
	if (in == NULL) {
		return NULL;
	}
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
	new_event->deposited	= source->deposited;
	new_event->location		= source->location;
	new_event->tof			= source->tof;
	new_event->particle		= source->particle;
	strncpy(new_event->orgin, source->orgin, ORIGIN_BUFFER);
	new_event->id			= source->id;
	return new_event;
}



// frees event malloc, returns NULL. For fmap
void* delete_event(void* in) {
	if (in == NULL) {
		return NULL;
	}
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
		// history_num = 0;
		previous_event = f(source);
		if (previous_event == NULL) {
			// probably at EOF, in any case we need to be done
			if (GENERAL_DEBUG) {
				printf("NULL event reached, ending\n");
			}
			return NULL;
		}
	}
	history_num = previous_event->number;
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
	return history;
}

// now for a quick array sorting implemetation

int scatter_dep_compare(void* va, void* vb) {
	if (va == NULL || vb == NULL) {
		return 0;
	}
	scatter* a = (scatter*)va;
	scatter* b = (scatter*)vb;
	if (a->deposit < b->deposit) {
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
	for (int j = low; j < high; j++) {
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
 * vec_to_phi
 * Takes a vector and returns the angle in phi as a double. Dead simple for use
 * by the trigger
 * Returns the angle in radians (from 0 to 2*pi), or an error value of -1
 */
double vec_to_phi(vec3d a) {
	if ((a.x == 0.0) && (a.y == 0.0)) {
		// cannot find phi as we are on the z axis, it is undefined. Return error
		return -1;
	}
	double angle = atan2(a.y, a.x);
	if (angle < 0) {
		angle += 2 * PI;
	}
	return angle;
}

int test_vec_to_phi() {
	vec3d case_a;
	vec3d case_b;
	vec3d case_c;
	int success = 0;
	case_a.x = 2.0;
	case_a.y = 0.0;
	double result = vec_to_phi(case_a);
	int run = double_equality(0.0, result, 0.0001);
	if (run != 1) {
		fprintf(stderr, "test_vec_to_phi: expected 0.0, got %lf\n",result);
	}
	success += run;
	case_b.x = 0.0;
	case_b.y = 1.5;
	result = vec_to_phi(case_b);
	run = double_equality(PI * 0.5, result, 0.0001);
	if (run != 1) {
		fprintf(stderr, "test_vec_to_phi: expected pi/2, got %lf\n",result);
	}
	success += run;
	case_c.x = 1.0;
	case_c.y = -1.0;
	result = vec_to_phi(case_c);
	run = double_equality(PI * 1.75, result, 0.0001);
	if (run != 1) {
		fprintf(stderr, "test_vec_to_phi: expected 5pi/4, got %lf\n",result);
	}
	success += run;
	if (success != 3) {
		fprintf(stderr, "test_vec_to_phi: tests failed.\n");
	}
	return success;
}

/*
 * phi_trigger:
 * Takes a scatter and the number of modules in phi that the detector has and
 * returns which module the scatter triggered. If it failed to trigger a module
 * 0 is returned. Otherwise, the number of the module as a binned value is
 * returned with 1 as the first module and the highest number of modules as the
 * final one.
 */
int phi_trigger(scatter* a, uint modules) {
	if (a == NULL) {
		return 0; // no trigger, but also lets not break anything
	}
	if (a->deposit < E_TRIGGER) {
		return 0;
	}
	double phi = vec_to_phi(a->loc);
	double segment = phi * (1.0 / (2 * PI)) * modules;
	int final = ((int)segment) + 1;
	return final;
}

// a simple function for calculating factorials for optimization check
uint factorial(uint a) {
	uint val = 1;
	for (; a > 0; a--) {
		val *= a;
	}
	return val;
}

// test factorial
int test_factorial() {
	uint result = factorial(5);
	if (result != 120) {
		fprintf(stderr, "test_factorial: failed, got %u\n",result);
		return 0;
	}
	return 1;
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
	double radius2 = event->location.x * event->location.x 
						+ event->location.y * event->location.y;
	// printf("particle %i, distance^2 %d", event->particle, event->count);
	if (sqrt(radius2) < distance) {
		if (fabs(event->location.z) < height){
			return 1;
		}
		return 0;
	}
	return 0;
}

// returns 1 for finding the first gamma it runs across has an IPS. Returns a 
// 2 for findign the second gamma it runs across has an IPS. These are then added
// to give 0,1,2,3 for the four possible cases: 0: no IPS, 1: only the first gamma,
// 2: only the second gamma, and 3: both gammas
int rec_in_paitent(llist* list, int reject_id) {
	if (list == NULL) {
		return 0;
	}
	// finds the first two gammas of the given history and looks to see if they
	// are just about 511 keV
	event* current = (event*)list->data;
	if ((current->particle == 22) && (reject_id == -1) && (current->id < 4)) {
		// we have our first gamma
		if (fabs(current->energy - ELECTRON_MASS) > 0.1) {
			// we had an IPS, so lets say so
			return 1 + rec_in_paitent(list->down, current->id);
		} else {
			reject_id = current->id;
			return rec_in_paitent(list->down, reject_id);
		}
	}
	if ((current->particle == 22) && (current->id != reject_id) && (current->id < 4)) {
		// different gamma!
		if (fabs(current->energy - ELECTRON_MASS) > 0.1) {
			// we had an IPS, so lets say so
			return 2;
		} else {
			return 0;
		}
	}
	return rec_in_paitent(list->down, reject_id);
}

/* 
 * Checks to see if there has been an in paitent scattering event in
 * the given history. An event is defined as a gamma showing up in the
 * list within the given radius of the z axis, with the previous instance
 * of the gamma having an energy greater by at least ENG_RNG. It finds
 * this information by calling rec_in_paitent, which checks over the
 * list for a gamma that fits the above conditions
 */
int in_patient(llist* list) {
	if (list == NULL)
		fprintf(stderr, "in_patient passed NULL pointer\n");
	list = list_head(list);
	return rec_in_paitent(list, -1);
}

// calculates the crystal number that given vector should be associated with
// this is a signed index number. Negative numbers mean -Z values. Index changes
// by 1 with rotations in phi, changes by the number of crystals in the layer
// radially inward for increasing R, changes by the number of crystals in a Z
// slice. This version assumes the same number of crystals in phi at all R values
int pos_to_crystal_index(vec3d pos, int crystals_in_phi, double crystal_r, double z_thickness) {
	// calculate the crystals in a z slice
	pos.z += 200.0;
	int z_silce_crystals = crystals_in_phi * ((int)ceil((double)SCANNER_THICKNESS / crystal_r));
	int z_slice_index = (int)floor(pos.z / z_thickness) * z_silce_crystals;
	// printf("pos_to: phi = %lf\n", vec_to_phi(pos));
	int phi_index = (int)floor((vec_to_phi(pos) / (PI * 2.0)) * crystals_in_phi);
	pos.z = 0.0;
	// printf("pos_to: r = %lf\n", vec_mag(pos));
	int r_index = (int)floor((vec_mag(pos) - (double)BORE_RADIUS) / crystal_r) * crystals_in_phi;
	return phi_index + z_slice_index + r_index;
}


vec3d crystal_index_to_pos(int index, int crystals_in_phi, double crystal_r, double z_thickness) {
	int z_silce_crystals = crystals_in_phi * ((int)ceil((double)SCANNER_THICKNESS / crystal_r));
	double z_dist = (index / z_silce_crystals) * z_thickness + (0.5 * z_thickness);
	int index_no_z = index % z_silce_crystals;
	double radius = (index_no_z / crystals_in_phi) * crystal_r + BORE_RADIUS + (0.5 * crystal_r);
	int index_no_r = index_no_z % crystals_in_phi;
	double phi = (((double)index_no_r + 0.5) / (double)crystals_in_phi) * (PI * 2.0);
	vec3d cartesian;
	cartesian.z = z_dist - 200.0;
	cartesian.x = radius * cos(phi);
	cartesian.y = radius * sin(phi);
	// printf("to_pos: phi = %lf\n", phi);
	// printf("to_pos: r = %lf\n", radius);
	return cartesian;
}

int test_crystal_indicies() {
	int pass = 0;
	// first test case, check a vector inside the detector and see if when converted back to
	// a vector it is close to the correct position
	int crystals_in_phi = 1024; // these values roughly match uEXPLORER
	double crystal_height = 1.5;
	double z_thickness = 0.276;
	vec3d initial_pos = three_vec(35.1, 34.3, -101.1);
	// this was just a weird position to choose.
	int index = pos_to_crystal_index(initial_pos, crystals_in_phi, crystal_height, z_thickness);
	vec3d final_pos = crystal_index_to_pos(index, crystals_in_phi, crystal_height, z_thickness);
	if (vec_mag(vec_sub(initial_pos, final_pos)) > 1.0) {
		// currently set to fail so that it can be checked by eye
		fprintf(stderr, "test_crystal_indicies-TEST FAILED: expected these two positions to be close:\n");
		fprintf(stderr, "\tinitial position: ");
		vec_print(initial_pos, stderr);
		fprintf(stderr, "\n\tend position: ");
		vec_print(final_pos, stderr);
		fprintf(stderr, "\n\tindex %i\n", index);
		pass += 1;
	}


	if (pass != 0) {
		fprintf(stderr, "!!! test_crystal_indicies: TESTS FAILED !!!\n\t%i tests failed\n", pass);
		return 1;
	}
	return 0;
}

void print_crystal(FILE* output, crystal* a) {
	fprintf(output, "E = %lf keV, t = %lf ns, index %i at ", a->energy, a->time, a->index);
	vec_print(crystal_index_to_pos(a->index, CRYSTAL_IN_PHI, CRYSTAL_HEIGHT, CRYSTAL_Z_THICKNESS), output);
	fprintf(output, ", with true N of");
	uint true_n_mask = 0b1;
	for (int i = 0; i < sizeof(int) * 8; i++) {
		if (true_n_mask & a->true_n) {
			fprintf(output, " %i,", i);
		}
		true_n_mask = true_n_mask << 1;
	}

	fprintf(output, " for gamma %i\n", a->true_gamma);
}

void* mappable_crystal_print(void* in) {
	print_crystal(stdout, (crystal*)in);
	return in;
}

// /*
//  * find_annihilation_point
//  * Takes a history and finds where the annihilation occured. If no
//  * annihilation occured it returns NULL. If one did occur it returns an
//  * array of doubles of length 3. This has the x,y,z position of the final
//  * location of the positron
//  * If it does not find an annihilation point it returns NULL
//  */
// vec3d* find_annihilation_point(llist *history) {
// 	if (history == NULL) {
// 		return NULL;
// 	}
// 	history = list_head(history);
// 	// looks for the first event that refers to a positron
// 	while (((event*)(history->data))->particle != -11) {
// 		history = history->down;
// 		if (history == NULL) {
// 			return NULL;
// 		}
// 	}
// 	// now look for the last event that refers to a positron
// 	while ((history->down != NULL) && (((event*)(history->down->data))->particle == -11)) {
// 		history = history->down;
// 	}
// 	// now at the last event, return the location as a 3-vector
// 	return vec_copy(((event*)(history->data))->location);
// }

/*
 * line_to_dot_dist:
 * finds the distance from a line defined by the start and end points to a
 * given point. The distance is how far the minimum distance is.
 */
double line_to_dot_dist(vec3d start, vec3d end, vec3d point) {
	vec3d num_first_term = vec_sub(start, end);
	vec3d num_sec_term = vec_sub(start, point);
	// vec3* denom_vec = vec_sub(end, start);
	double numerator = vec_mag(vec_cross(num_first_term, num_sec_term));
	double denomenator = vec_mag(num_first_term);
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
	double b_der = .5 * (first_term + 1.);
	uncert_b = uncert_b * b_der;
	if (GENERAL_DEBUG) {
		printf("expected_uncert_b: b_der = %lf\n", b_der);
		printf("expected_uncert_b: first_term = %lf\n", first_term);
		printf("expected_uncert_b: energy uncertainty is %lf\n", uncert_b);
		printf("expected_uncert_b: theta uncertainty is %lf\n", uncert_theta);
	}

	// double combined = 
	// if (combined < 0) {
	// 	fprintf(stderr, "expected_uncert_b: combined negative.\n\tHow did you get here?\n");
	// 	fprintf(stderr, "\tcombined = %lf\n\tuncert_theta = %lf\n", combined, uncert_theta);
	// 	fprintf(stderr, "\tuncert_b = %lf\n", uncert_b);
	// }
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
	vec3d ab = vec_sub(b->loc, a->loc);
	// get the vector from b->c
	vec3d bc = vec_sub(c->loc, b->loc);
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
		if (delta_e < 0) {
			fprintf(stderr, "expected_energy_b: uncertanty negative: %lf\n",delta_e);
		}
		uncert[0] = delta_e;
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

	vec3d point_a = three_vec(0., 0., 0.);
	vec3d point_b = three_vec(0., 3., 0);
	vec3d point_c = three_vec(0., 2., 1.73205);
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

// dot product of vector first->second locations and the electron direction
// at second
double scatter_dir_dot(scatter* first, scatter* second) {
	vec3d a_less_b = vec_sub(first->loc, second->loc);
	double dot = vec_dot(a_less_b, second->dir);
	return dot;
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


/*
 * closest_gamma:
 * finds the id of the closest gamma to the given location. The gamma is sourced
 * from a given history. The gamma is required to be from an annihilation
 * If the search fails returns zero. We will only check for gammas labeled 2 and 3
 * This will hurt the search sometimes but will cutout wierd PE things
 */
int closest_gamma(llist* history, vec3d target, double* deposited) {
	if (history == NULL) {
		return 0;
	}
	history = list_head(history);
	double distance = INFINITY;
	int best_id = 0;
	double deposit = 0;
	event* working_event = NULL;
	while (history != NULL) {
		// first, define our working event
		working_event = (event*)history->data;
		if ((working_event->particle == 22)) {
			// it is a gamma from an annihilation
			double check_dist = vec_dist(target, working_event->location);
			if ((distance > check_dist) && (working_event->id < 4)) {
				distance = check_dist;
				best_id = working_event->id;
				deposit = working_event->deposited;
			}
		}
		// move on to next event
		history = history->down;
	}
	if (deposited != NULL) {
		deposited[0] = deposit;
	}
	return best_id;
}

/*
 * build_scatters
 * takes the list of interactions in the detector volume and makes a list of
 * energy deposition locations that are the scatter locations. This is the
 * current version of the "cluster finding" algorithm. A list of all of the
 * scatters is passed out, not combined into crystals. No resolutions are
 * applied here.
 */
llist* build_scatters(llist* detector_history) {
	if (detector_history == NULL) {
		return NULL;
	}
	// looks to find the first instance of each electron showing up. Using the
	// location of that a scatter
	// of a gamma with the same identifier as id. If all is good the location is
	// added as a scatter with the electron's energy plus any energy deposited by
	// a gamma at that location

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
						// the current electron is in our list of handled electrons
						electron_checked = 1;
					}
					check = check->down;
				}
			}
			// electron_checked now contains if we have previously looked at this
			// electron
			if (!electron_checked) {
				//  so which gamma is it from?
				double deposit;
				int gamma_id = closest_gamma(detector_history, ((event*)detector_history->data)->location, &deposit);
				// time to add this electron to the scatter list

				vec3d scatter_loc = (((event*)detector_history->data)->location);
				event* cur_event = (event*)detector_history->data;

				// char has_dir = 0; // we aren't doing anything with electrons so let's not pretend
				scatter* add_scatter;
				double photons = cur_event->energy * P_PER_KEV;
				add_scatter = new_scatter(scatter_loc, three_vec(0,0,0), gamma_id, cur_event->energy + deposit, cur_event->tof, sqrt(photons) / P_PER_KEV, SPC_UNCERT_PLANE, TIME_UNCERT_CM);
				scatter_truth* truth_info = (scatter_truth*)malloc(sizeof(scatter_truth));
				truth_info->true_eng = cur_event->energy;
				truth_info->true_time = cur_event->tof;
				add_scatter->truth = truth_info;

				scatter_list = add_to_bottom(scatter_list, add_scatter);

				// add the electron to the list of electrons we have checked
				int* electron_id = (int*)malloc(sizeof(int));
				*electron_id = ((event*)detector_history->data)->id;
				checked_electrons = add_to_bottom(checked_electrons, electron_id);
			}
		}
		detector_history = detector_history->down;
	}

	fmap(checked_electrons, free_null);
	delete_list(checked_electrons);
	// go through and check which scatter N each one is
	int gamma_2_N = 0;
	int gamma_3_N = 0;
	for (llist* ordering = list_tail(scatter_list); ordering != NULL; ordering = ordering->up) {
		scatter* working_scatter = (scatter*)(ordering->data);
		if ((working_scatter->truth != NULL) && (working_scatter->has_dir == 2)) {
			gamma_2_N += 1;
			working_scatter->truth->true_n = gamma_2_N;
		} else if ((working_scatter->truth != NULL)) {
			gamma_3_N += 1;
			working_scatter->truth->true_n = gamma_3_N;
		}
	}

	// return the list of scatters for consolidation into specific crystals

	return scatter_list;
}

// takes a list of scatters and makes a list of crystals that lit up, and how
// much they did so.
llist* build_crystals(llist* scatters) {
	if (scatters == NULL) {
		return NULL;
	}
	llist* list_of_crystals = NULL;

	// we will work our way down the list of scatters adding each scatter's
	// energy to the right crystal. If a scatter has the same crystal ID as
	// an already existing scatter then we do not make a new crystal and instead
	// add the energy to the current one
	scatters = list_head(scatters);
	while (scatters != NULL) {
		scatter* working_scatter = (scatter*)(scatters->data);
		int scatter_index = pos_to_crystal_index(working_scatter->loc, CRYSTAL_IN_PHI, CRYSTAL_HEIGHT, CRYSTAL_Z_THICKNESS);

		// check to see if the crystals list has the current index
		llist* check_crystals = list_head(list_of_crystals);
		int preexisting_crystal = 0;
		while ((check_crystals != NULL) && (!preexisting_crystal)) {
			crystal* working_crystal = (crystal*)(check_crystals->data);
			if (working_crystal->index == scatter_index) {
				working_crystal->energy += working_scatter->deposit;
				if (working_crystal->time > working_scatter->time) {
					working_crystal->time = working_scatter->time;
				}
				if (working_scatter->truth != NULL) {
					working_crystal->true_n = working_crystal->true_n | ( 0b1 << working_scatter->truth->true_n);
				}
				preexisting_crystal = 1;
			}


			check_crystals = check_crystals->down;
		}
		if (preexisting_crystal == 0) {
			// no preexisting crystal to add the energy to, need to make one;
			crystal* new_crystal = (crystal*)malloc(sizeof(crystal));
			new_crystal->energy = working_scatter->deposit;
			new_crystal->index = scatter_index;
			new_crystal->pos = crystal_index_to_pos(scatter_index, CRYSTAL_IN_PHI, CRYSTAL_HEIGHT, CRYSTAL_Z_THICKNESS);
			new_crystal->time = working_scatter->time;
			new_crystal->true_gamma = working_scatter->has_dir;
			if (working_scatter->truth != NULL) {
				new_crystal->true_n = 0b1 << working_scatter->truth->true_n;
			} else {
				new_crystal->true_n = 0b0;
			}
			list_of_crystals = add_to_top(list_of_crystals, new_crystal);
		}

		scatters = scatters->down;
	}

	// here is where we would add energy and time resolution randomness.
	for (llist* unrandom = list_head(list_of_crystals); unrandom != NULL; unrandom = unrandom->down) {
		crystal* working_crystal = (crystal*)(unrandom->data);

		double time_var = 0.0;
		double eng_var = 0.0;
		for (int i = 0; i < UNCERT_REP; i++) {
			// creates a randomly distributed value
			time_var += drand48();
			eng_var += drand48();
		}
		time_var -= ((float)(UNCERT_REP) * 0.5);
		eng_var  -= ((float)(UNCERT_REP) * 0.5);
		// distance and time now have their variation size

		// calculate the amount to adjust the time by
		double time_adjust = time_var * (TIME_UNCERT_CM / SPD_LGHT);
		// blur the time
		working_crystal->time += time_adjust;

		// calculate the new energy using the number of photons we will count
		double new_eng = (eng_var * (sqrt(working_crystal->energy * P_PER_KEV) / P_PER_KEV)) + working_crystal->energy;
		working_crystal->energy = new_eng; // set the new energy

	}

	return list_of_crystals;
}

// takes a crystal and converts it to a scatter for use by a search algorithm
scatter* crystal_to_scatter(crystal* in) {
	if (in == NULL) {
		return NULL;
	}
	scatter* new = new_scatter(in->pos, three_vec(0,0,0), in->true_gamma, in->energy, in->time, sqrt(in->energy / P_PER_KEV), SPC_UNCERT_PLANE, TIME_UNCERT_CM);
	new->truth = (scatter_truth*)malloc(sizeof(scatter_truth));
	int true_n_mask = 0b1;
	for (int i = 0; i < sizeof(int) * 8; i++) {
		if (true_n_mask & in->true_n) {
			new->truth->true_n = i;
			break;
		}
		true_n_mask = true_n_mask << 1;
	}
	// new->truth->true_n = in->true_n;
	return new;
}

// takes a list of crystal interactions and makes a pair of lists of scatters for processing
// in the standard inverse_kinematics methods.
// The two lists are made by cluster finding to determine the sets of spatially close scatters.
// The cluster finding is done by k-means. Choose two points as inital means. Then for each crystal
// calculate the closer mean and move the means to the average location of their crystals

// if the code fails for some reason the function returns a 0. Otherwise it returns a 1.
int cluster_finding(llist* list_of_crystals, llist** cluster_1, llist** cluster_2) {
	if ((list_of_crystals == NULL) || (cluster_1 == NULL) || (cluster_2 == NULL)) {
		// not all the arguments are available, cannot continue
		return 0;
	}
	int length = list_length(list_of_crystals);
	if (length < 2) {
		return 0;
	} else if (length == 2) {
		// make one scatter list from each crystal
		scatter* first = crystal_to_scatter((crystal*)(list_of_crystals->data));
		cluster_1[0] = add_to_top(NULL, first);
		scatter* second = crystal_to_scatter((crystal*)(list_of_crystals->down->data));
		cluster_2[0] = add_to_top(NULL, second);
		return 1;
	}

	// OK we now need to do actual cluster finding
	// set the two initial mean locations to be the first and last crystals in the list
	// for prior data handling reasons these are likely to be close to the correct positions
	vec3d mean_1 = ((crystal*)(list_head(list_of_crystals)->data))->pos;
	vec3d mean_2 = ((crystal*)(list_tail(list_of_crystals)->data))->pos;

	// do Lloyd's algorithm for x iterations

	for (int i = 0; i < 10; i++) {
		llist* set_of_crystals = list_head(list_of_crystals);
		vec3d adj_mean_1 = three_vec(0,0,0);
		int count_mean_1 = 0;
		vec3d adj_mean_2 = three_vec(0,0,0);
		int count_mean_2 = 0;
		while (set_of_crystals != NULL) {
			crystal* working_crystal = (crystal*)(set_of_crystals->data);
			// check to see which mean is closer
			vec3d from_mean_1 = vec_sub(working_crystal->pos, mean_1);
			vec3d from_mean_2 = vec_sub(working_crystal->pos, mean_2);
			if (vec_mag(from_mean_1) < vec_mag(from_mean_2)) {
				// mean_1 was closer
				count_mean_1 += 1;
				adj_mean_1 = vec_add(adj_mean_1, from_mean_1);
			} else {
				// mean_2 was closer
				count_mean_2 += 1;
				adj_mean_2 = vec_add(adj_mean_2, from_mean_2);
			}

			set_of_crystals = set_of_crystals->down;
		}
		// we have now found what cluster each is in and how much to move
		// the means by
		adj_mean_1 = vec_scaler(adj_mean_1, 1.0 / (double)count_mean_1);
		adj_mean_2 = vec_scaler(adj_mean_2, 1.0 / (double)count_mean_2);
		mean_1 = vec_add(mean_1, adj_mean_1);
		mean_2 = vec_add(mean_2, adj_mean_2);
	}

	// we now have our set mean positions, lets make the two scatter lists based on that
	llist* set_of_crystals = list_head(list_of_crystals);
	cluster_1[0] = NULL;
	cluster_2[0] = NULL;
	// set up some truth checking data to know how well the cluster finding went
	int correct_find = 1;
	int first_gamma = 0;
	int second_gamma = 0;

	while (set_of_crystals != NULL) {
		crystal* working_crystal = (crystal*)(set_of_crystals->data);
		scatter* new = crystal_to_scatter(working_crystal);
		// check which cluster this one is in
		if (vec_mag(vec_sub(working_crystal->pos, mean_1)) < vec_mag(vec_sub(working_crystal->pos, mean_2))) {
			// this was placed into cluster 1
			cluster_1[0] = add_to_bottom(cluster_1[0], new);
			if (first_gamma == 0) {
				first_gamma = working_crystal->true_gamma;
			} else {
				if (first_gamma != working_crystal->true_gamma) {
					correct_find = 0;
				}
			}
		} else {
			cluster_2[0] = add_to_bottom(cluster_2[0], new);
			if (second_gamma == 0) {
				second_gamma = working_crystal->true_gamma;
			} else {
				if (second_gamma != working_crystal->true_gamma) {
					correct_find = 0;
				}
			}
		}

		set_of_crystals = set_of_crystals->down;
	}
	if (correct_find) {
		correct_clustering += 1;
	}
	total_clustering += 1;
	// we have now divided the two sets of scatters into clusters

	return 1;
}


/* find the electron paths that were recorded. If there are any number of
 * electrons than two then exit returning NULL. If only two electrons are there
 * then make them into scatters, doing all the good stuff with adding randomness
 * This randomness is not even in all directions. It has the FWHM in the plane
 * and radial directions as given by the constants at the top of the file.
 * If both energies after randomness in energy is applied pass the energy cut
 * then the two locations are passed out as endpoints with whatever timing
 * resolution
 */

scatter** find_endpoints(llist* detector_history, double energy_percent) {
    if (detector_history == NULL) {
        return NULL;
    }
    detector_history = list_head(detector_history);
	llist* search_history = detector_history;
	int first_id = 0;
	int second_id = 0;

	while (search_history != NULL) {
		// find the lowest two ID gammas (they will be from the annihilation)
		event* current_event = (event*)search_history->data;
		if ((current_event->particle == 22) && (current_event->id != first_id) && (current_event->id != second_id)) {
			// we have a gamma, and it isn't the same as the (possibly alreadly found)
			// first gamma
			if (first_id == 0) {
				first_id = current_event->id;
			} else {
				if ((current_event->id < second_id) || (second_id == 0)) {
					second_id = current_event->id;
				}
			}
		}
		search_history = search_history->down;
	}
	if (second_id == 0) {
		// no second gamma was found, so no line of response can be made.
		return NULL;
	}
	// printf("gamma id 1: %i\tgamma id 2: %i\n", first_id, second_id);

	llist* scat_list1 = NULL;
	llist* scat_list2 = NULL;
	llist* full_scat_list = build_scatters(detector_history);
	llist* crystal_list = build_crystals(full_scat_list);
	if (CLUSTER_DEBUG) {
		printf("crystal list: \n");
		fmap(crystal_list, mappable_crystal_print);
	}
	int success = cluster_finding(crystal_list, &scat_list1, &scat_list2);
	if (CLUSTER_DEBUG) {
		printf("scatter list 1: \n");
		fmap(scat_list1, mappable_print_scatter);
		printf("scatter list 2: \n");
		fmap(scat_list2, mappable_print_scatter);
	}

	fmap(full_scat_list, delete_scatter);
	delete_list(full_scat_list);
	fmap(crystal_list, free_null);
	delete_list(crystal_list);
	if (!success) {
		if (scat_list1 != NULL) {
			fmap(scat_list1, delete_scatter);
			delete_list(scat_list1);
		}
		if (scat_list2 != NULL) {
			fmap(scat_list2, delete_scatter);
			delete_list(scat_list2);
		}
		return NULL;
	}

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

    // we now have a list of scatter locations. Check if there is a scatter that passes the energy cut for the first one
	scatter* first_endpoint = NULL;
	llist* current_location = list_head(scat_list1);
	while ((current_location != NULL) && (first_endpoint == NULL)) {
		// printf("first scatter list energy error: %lf\n", fabs(((scatter*)(current_location->data))->deposit - ELECTRON_MASS) / ELECTRON_MASS);
		// printf("first scatter list energy: %lf\n", ((scatter*)(current_location->data))->deposit);
		if (fabs(((scatter*)(current_location->data))->deposit - ELECTRON_MASS) / ELECTRON_MASS < energy_percent) {
			first_endpoint = (scatter*)(current_location->data);
		} else {
			current_location = current_location->down;
		}
	}
	scatter* second_endpoint = NULL;
	current_location = list_head(scat_list2);
	while ((current_location != NULL) && (second_endpoint == NULL)) {
		// printf("second scatter list energy error: %lf\n", fabs(((scatter*)(current_location->data))->deposit - ELECTRON_MASS) / ELECTRON_MASS);
		// printf("second scatter list energy: %lf\n", ((scatter*)(current_location->data))->deposit);
		if (fabs(((scatter*)(current_location->data))->deposit - ELECTRON_MASS) / ELECTRON_MASS < energy_percent) {
			second_endpoint = (scatter*)(current_location->data);
		} else {
			current_location = current_location->down;
		}
	}
	if ((first_endpoint != NULL) && (second_endpoint != NULL)) {
		
		scatter** return_vals = (scatter**)malloc(2 * sizeof(scatter*));
		return_vals[0] = copy_scatter(first_endpoint);;
		return_vals[1] = copy_scatter(second_endpoint);

		fmap(scat_list1, delete_scatter);
		fmap(scat_list2, delete_scatter);
		delete_list(scat_list1);
		delete_list(scat_list2);

		return return_vals;
	}


    fmap(scat_list1, delete_scatter);
    fmap(scat_list2, delete_scatter);
    delete_list(scat_list1);
    delete_list(scat_list2);
    return NULL;
}



// /* 
//  * first_scat_miss_transverse:
//  * Takes the two endpoints and a 3-vector defining the location of the
//  * annihilation and finds the distance between the line-of-responce and the
//  * annihilation. If there is a problem (such as not having two endpoints) it
//  * returns -1 as a reject value.
//  */
// double first_scat_miss_transverse(scatter** endpoints, vec3d* annh_loc) {
// 	if ((endpoints == NULL) || (annh_loc == NULL)) {
// 		return -1;
// 	}
// 	if ((endpoints[0] == NULL) || (endpoints[1] == NULL)) {
// 		return -1;
// 	}
// 	// find the scattering point of the first scatter

// 	vec3d* first_spot = vec_copy(endpoints[0]->loc);

// 	vec3d* second_spot = vec_copy(endpoints[1]->loc);

// 	return line_to_dot_dist(first_spot, second_spot, annh_loc);
// }

/* first_scat_miss_transverse:
 * takes a LOR and the location of the annihilation and returns the transverse
 * distance between the two. That is: the distance between the center of the LOR
 * and the annihilation location perpendicular to the direction of the LOR
 */
double first_scat_miss_transverse(lor* lor, vec3d annh_loc) {
	if (lor == NULL) {
		return -1;
	}
	vec3d offset = vec_sub(annh_loc, lor->center);
	vec3d perpendicular = vec_cross(offset, lor->dir);
	double dist = vec_mag(perpendicular);
	return dist;
}

/* first_scat_miss_longitudinal:
 * takes a LOR and the location of the annihilation and returns the
 * longitudinal distance between the two. That is: the distance between the
 * center of the LOR and the annihilation projected onto the direction of the
 * LOR.
 */
double first_scat_miss_longitudinal(lor* lor, vec3d annh_loc) {
	if (lor == NULL) {
		return -1;
	}
	vec3d offset = vec_sub(annh_loc, lor->center);
	double offset_dist = vec_mag(offset);
	if (offset_dist < ENG_RNG) {
		return offset_dist;
	}
	double dist = vec_dot(offset, lor->dir);
	return dist;
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
	if (item->particle == 22) {
		return item->id;
	}
	
	// no longer works without origin data
	// if (item->orgin[0] == 'a') {
	// 	return item->id;
	// }
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
	vec3d center_subtraction = vec_sub(a->loc, b->loc);
	vec3d center_half = vec_scaler(center_subtraction, 0.5);
	vec3d geometric_center = vec_add(b->loc, center_half);
	vec3d ba_unit = vec_norm(center_half);
	double time_delta = 0.5 * (b->time - a->time);
	vec3d displacement = vec_scaler(ba_unit, SPD_LGHT * time_delta);
	lor* new = (lor*)malloc(sizeof(lor));
	new->center = vec_add(geometric_center, displacement);
	new->dir = ba_unit;
	double a_space = a->space_uncert;
	double b_space = b->space_uncert;
	if (a_space < 0) {
		a_space = .1;
	}
	if (b_space < 0) {
		b_space = 0.1;
	}

	new->transverse_uncert = (CRYSTAL_Z_THICKNESS) * 0.5;
	new->long_uncert = sqrt(a->time_uncert * a->time_uncert + b->time_uncert * b->time_uncert) + (a_space + b_space) / 2;

	return new;
}

void print_lor(FILE* output, lor* lor) {
	fprintf(output, "%f, %f, %f,", lor->center.x, lor->center.y, lor->center.z);
	fprintf(output,  " %f, %f, %f,", lor->dir.x, lor->dir.y, lor->dir.z);
	fprintf(output, " %f, %f", lor->long_uncert, lor->transverse_uncert);
}


int main(int argc, char **argv) {
	// quick reminder that argc is the number of arguments,
	// argv points to the strings of the arguments

	// we are looking for an input histories file, an output file name,
	// and an inner radius of the detector this may later be changed to
	// input histories, input data file, output file name
	if (argc < 4) {
		printf("Unable to run. Expected at least 3 arguments got %i.\n", argc - 1);
		printf("Expects an input detector volume ");
		printf("history file, output file name, and");
		printf(" energy cutoff percent. Can be followed by -b for reading binary files\n");
		return 1;
	}
	// FILE* in_histories = fopen(argv[1], "r");
	// if (in_histories == NULL) {
	// 	printf("Unable to open history file\n");
	// 	return 1;
	// }
	FILE* in_det_histories;
	int binary = 0;
	if (argc > 4) {
		for (int i = 1; i < argc; i++) {
			if (!strcasecmp(argv[i], "-b")) {
				// binary flag
				binary = 1;
			}
		}
	}
	// if (binary) {
	in_det_histories = fopen(argv[1], "r");
	// } else {
		// in_det_histories = fopen(argv[1], "r");
	// }

	// FILE* in_det_histories = fopen(argv[1], "r");
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

	double energy_cutoff = strtod(argv[3], NULL);

	if ((!test_expected_energy()) || (test_vec_to_phi() != 3)) {
		fprintf(stderr, "tests failed, exiting\n");
		return 1;
	}
	llist *in_det_hist = NULL;
	if (binary) {
		in_det_hist = load_historyb(in_det_histories, read_line_binary);
	} else {
		in_det_hist = load_historyb(in_det_histories, read_line);
	}
	uint run_num = 0;
	first_scat_hypot = 0;
	second_scat_hypot = 0;

	// if (pthread_mutex_init(&file_lock, NULL)) {
	// 	fprintf(stderr, "pthread: unable to make volume lock, exiting.\n");
	// 	return(1);
	// }

	// for (int i = 0; i < MAX_THREAD_CALLS; i++) {
	// 	// set the array of threads to be empty
	// 	tid[i] = -1;
	// }
	// int cur_thread = 0;


	// begin the primary loop over all histories
	while (in_det_hist != NULL) {
		// print an update to how far we have made it
		if ((run_num / 1000000) * 1000000 == run_num) {
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

		scatter** endpoints = find_endpoints(in_det_hist, energy_cutoff);

		if (endpoints == NULL) {

		} else {
			if (in_patient(in_det_hist)) {
				scatter_counts++;
			} else {
				true_counts++;
			}
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
		if (binary) {
			in_det_hist = load_historyb(in_det_histories, read_line_binary);
		} else {
			in_det_hist = load_historyb(in_det_histories, read_line);
		}		run_num++;

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

	printf("Scattered: %u\n", scatter_counts);
	printf("True: %u\n", true_counts);

	printf("Cluster finding performace: %lf %%\n", 100. * ((float)correct_clustering / (float)total_clustering));


	fclose(in_det_histories);
	fclose(lor_output);
	return 0;
}

