#include "truth_assign.h"
#include <stdio.h>
#include <string.h>
#include "llist.h"
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"
#include "ray_trace.h"

#define ENG_RNG 0.001
#define COMP_INT 1667457891
#define ELECTRON_MASS 510.999
#define SPD_LGHT 29.98
#define PI 3.141592653589
#define FIRST_N 5
#define LARGEST 10 // max scatters in a tree
#define SKIP 0 // number of scatters to skip at end of tree
#define KEEP_SINGLES 1
#define MAX_SINGLE_FOM 3.0 // max sigma for a scatter to be accepted as possible
#define MIN_SCAT_ENG 0.0 //10.0 // minumum energy for a scatter to be observed
#define E_TRIGGER 0.0 //20.0 // energy in keV to hit trigger 
#define PHI_MODULES 12 // modules for trigger
#define MODULE_SEPERATION 3 // min number of modules separation for trigger
#define NEVER_CUT 0
#define BINDING_WIDTH 5.0 // a 1 sigma width in keV for the effects of electron binding energy and doppler
#define USE_SCATTER_DISTANCE 0 // flag on if to use the distance between scatters as part of FOM
#define BORE_RADIUS 45 // measured in cm
#define SCANNER_THICKNESS 15 // measured in cm
#define CRYSTAL_IN_PHI 1024 // these values roughly match uEXPLORER
#define CRYSTAL_HEIGHT 1.5
#define CRYSTAL_Z_THICKNESS 0.276


double time_uncert_cm = 6.36; // in cm for one sigma, NOT ps or ns FWHM
double spc_uncert = 0.1; // cm
double P_per_keV = 30.0; // photons/keV
int run_time_rand = 1;
#define UNCERT_REP 12
#define SPC_UNCERT_PLANE 0.3 // uncertainty in the plane of the bore
#define SPC_UNCERT_RAD 1.0  // uncertainty radial to the bore

#define READ_DEBUG 0
#define GENERAL_DEBUG 0
#define TREE_DEBUG 0
#define SCATTER_LIST_DEBUG 0
#define GRAPHVIZ_DEBUG 0
#define CLUSTER_DEBUG 0
FILE* debug_graphs = NULL;
FILE* debug_scatter_lists = NULL;
uint graph_id;
#define LOR_GROUP 0 // iff 1 outputs true and misID
#define CUT_IPS 0 // do not output to true/misID files if an IPS happened

//OUTDATED if 0 output all found lors, 1 only with both T=R=1, 2 
// only R!=T=1 for only 1, 3 only R!=T=1 for both 

int alpha_n_1[FIRST_N];
int alpha_n_2[FIRST_N];
float alpha_eng_distro_n_1[FIRST_N];
float alpha_eng_distro_n_2[FIRST_N];
float n_eng_distro_n_1[FIRST_N];
float n_eng_distro_n_2[FIRST_N];
uint hist_num;
uint total_scatters = 0;
uint path_scatters = 0;
unsigned long long possible_branches = 0;
uint missed_reconstructions = 0; // add 1 for each pair of gammas that made it to
// the detector but did not get reconstructed
uint missed_reconstruction_IPS_mask = 1; // stop the addition of a missed reconstruction
// if the reconstruction would have been of an IPS

int branch_scatters_1 = 0;
int branch_scatters_2 = 0;
float first_FOM = -1;
float second_FOM = -1;
int used_scatters_1 = -1;
int used_scatters_2 = -1;

int correct_clustering = 0;
int total_clustering = 0;

float bore_position[3] = {0,0,0};
float bore_dimentions[3] = {45, 200, 0};
shape* detector_bore = NULL;

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
	// previous_event = duplicate_event(previous_event);
	// WHOOPS THAT WAS A MEMORY LEAK!
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
	double angle = atan2(a.y,a.x);
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

// as copied from stack overflow
float myErfInv2(float x){
   float tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = logf(x);

   tt1 = 2/(PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
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

// /* 
//  * Checks to see if there has been an in paitent scattering event in
//  * the given history. An event is defined as a gamma showing up in the
//  * list within the given radius of the z axis, with the previous instance
//  * of the gamma having an energy greater by at least ENG_RNG. It finds
//  * this information by calling rec_in_paitent, which checks over the
//  * list for a gamma that fits the above conditions
//  */
// int in_patient(llist* list, double radius, double height, int id) {
// 	if (list == NULL)
// 		fprintf(stderr, "in_patient passed NULL pointer\n");
// 	list = list_head(list);
// 	return rec_in_paitent(list, radius, height, id);
// }


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
		// add in a binding energy uncertainty:
		delta_e = add_quadrature(delta_e, BINDING_WIDTH);
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
 * dist_probability:
 * Takes a pair of scatters from a detector and the predicted energy of the gamma
 * that went from one to the other. It then calculates the distance between the
 * scatters and using an approximation of the cross-section values for 
 */
double dist_sigma(scatter* a, scatter* b, double energy) {
	if ((a == NULL) || (b == NULL)) {
		return -INFINITY; // very bad, and noticably meaningless result (it's negative)
	}
	// To handle a gamma that passes back through the bore of the detector, subtract the distance
	// traveled through a cylinder with the inner radius of the detector.
	// To determine the distance travled within the detector medium calculate the distance between
	// the two scatters and then subtract the distance that ray travels inside the radius of the bore
	vec3d outgoing = vec_sub(b->loc, a->loc);
	double total_distance = vec_mag(outgoing);
	double crossing_distance = total_distance;
	// make a ray from loc to remaining
	// check the direction that the ray would have. If it is not inward we do not need to raytrace
	// dot the vector pointing from the center with the direction vector
	if (vec_dot(three_vec(outgoing.x, outgoing.y, 0), three_vec(a->loc.x, a->loc.y, 0)) < 0.0) {
		vec3d outgoing_unit = vec_scaler(outgoing, 1.0/total_distance);
		ray* outgoing_ray = ray_build(a->loc, outgoing_unit);
		traversal* bore_distance = exit_cyl(outgoing_ray, detector_bore, NULL);
		if (bore_distance != NULL) {
			crossing_distance += bore_distance->t;
			traversal_free(bore_distance);
		}
		ray_free(outgoing_ray);
	}
	// whats the probability this happened?
	// requires the interaction cross section of the material at various energies
	double cross_sections[10] = {0.166, 0.142, 0.128, 0.118, 0.109, 0.102, 0.096, 0.092, 0.088, 0.083};
	double energy_cross_section;
	if (energy < 50.0) {
		energy_cross_section = cross_sections[0];
	} else if (energy > 500.0) {
		energy_cross_section = cross_sections[9];
	} else {
		energy_cross_section = cross_sections[(int)round(energy / 50.0) - 1];
	}
	double probability = exp(-fabs(crossing_distance) * energy_cross_section); // distance / (1/e distance)
	// printf("distance: %lf,\t cross section: %lf, \t exponent = %lf\n",crossing_distance, energy_cross_section, -fabs(crossing_distance) * energy_cross_section);
	double sigma = myErfInv2(1.0 - probability);
	// printf("Prob: %lf, \t sigma %lf\n", probability, sigma);
	return sigma;
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
	if ((source == NULL) || (source_len <= 1)) {
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

double recursive_search(double best, double current, double inc_eng, double inc_uncert, scatter* origin, scatter* loc, scatter** remaining, uint remain_count, llist** path) {
	if ((loc == NULL) || (remaining == NULL) || (origin == NULL) || (remain_count <= SKIP)) {
		if (TREE_DEBUG) {
			printf("recursive_search: end of branch\n");
		}
		if (GRAPHVIZ_DEBUG) {
			fprintf(debug_graphs, "%u;\n", graph_id);
			fprintf(debug_graphs, "%u [label=\"end of check\"];\n", graph_id);
			graph_id++;
		}
		// if ((loc != NULL) && (path != NULL) && (loc->truth != NULL)) {
		// 	float* best_N_return = (float*)malloc(sizeof(float) * 2);
		// 	best_N_return[0] = loc->truth->true_n;
		// 	best_N_return[1] = loc->deposit;
		// 	path[0] = add_to_top(NULL, best_N_return);
		// }
		return current;
	}
	if (best < current) {
		// the best find is better than our current find, we are not doing well
		if (TREE_DEBUG) {
			printf("recursive_search: best better than current\n");
		}
		return current;
	}

	double OUT_OF_PLANE_WEIGHT = 0.; // how many sigma a 90 deg out of plane
	// scatter counts for 
	if (TREE_DEBUG) {
		printf("recursive_search: Origin N = ");
		if (origin->truth != NULL) {
			printf("%i, N = ", origin->truth->true_n);
		} else {
			printf("MISSING, N = ");
		}
		if (loc->truth != NULL) {
			printf("%i, ", loc->truth->true_n);
		} else {
			printf("MISSING, ");
		}
		printf("best = %lf, current = %lf\n", best, current);
	}
	uint graph_label;
	if (GRAPHVIZ_DEBUG && (loc->truth != NULL)) {
		fprintf(debug_graphs, "%u;\n", graph_id);
		fprintf(debug_graphs, "%u [label=\"N=%i\\nbest=%lf\\ncur=%lf\\neng=%lf keV\"];\n",graph_id, loc->truth->true_n,best,current, loc->deposit);
		graph_label = graph_id;
		graph_id++;
	}

	double energy_uncert;
	vec3d in = vec_sub(origin->loc, loc->loc);

	double better_find = INFINITY;

	// begin N vs alpha work
	llist* continuing_path = NULL;
	llist* best_path = NULL;
	// end N vs alpha work

	for (int i = 0; i < remain_count; i++) {
		// first calcuate the probability of our current deviation
		double prediction = expected_energy_b(origin, loc, remaining[i], &energy_uncert);
		if (energy_uncert < 0) {
			fprintf(stderr, "recursive_search: negative energy error found (%lf). Check given scatter errors\n",energy_uncert);
			fprintf(stderr, "\terror thrown on history %i\n", hist_num);
		}
		// calculate the combined error of the incoming energy and expectation
		double comb_uncert = add_quadrature(inc_uncert, energy_uncert);
		double energy_error = fabs(inc_eng - prediction) / comb_uncert; // error between
		// double new_energy = inc_eng - loc->deposit; // could maybe be the value
		// of expectation - loc->deposit
		double new_energy = inc_eng - loc->deposit;
		double new_eng_uncert = add_quadrature(inc_uncert, loc->eng_uncert);
		// the a->b gamma prediction and the value from earlier in the chain
		double step_error = energy_error;
		if ((step_error + current < best) && (step_error < MAX_SINGLE_FOM) && USE_SCATTER_DISTANCE) {
			// no sense in checking other constraints if we already fail
			// find the error based on the distance traveled between the two interaction locations
			double dist_error = dist_sigma(loc, remaining[i], prediction);
			step_error = add_quadrature(step_error, dist_error);
		}

		// now find the error in electron direction plane
		if (loc->has_dir) {
			// only do this if there is a point. No electron direction then you
			// can't check if it is in the plane.

			vec3d outgoing = vec_sub(loc->loc, remaining[i]->loc);
			// direction heading out from the current location to the next location
			vec3d scatter_plane_normal = vec_cross(in, outgoing);
			// creates a normal to the plane of the scattering origin->loc->i
			vec3d plane_norm_hat = vec_norm(scatter_plane_normal);
			// normalizes the plane perpendicular
			double out_of_plane = vec_dot(plane_norm_hat, loc->dir);
			// has the value of cosine of the angle between the plane normal and the
			// direction of the electron. This means that for perfect in the plane
			// behavior the value is 0, out of the plane +- 1.
			step_error = add_quadrature(step_error, abs(out_of_plane * OUT_OF_PLANE_WEIGHT));
		}
		// done calcuating what the current paths uncertanty is, lets combine it
		// with the current uncertanty to check if we are worse than the best
		// solution so far.
		// double combined_error = add_quadrature(step_error, current);
		double combined_error = step_error + current;
		uint part_graph_id;
		if(GRAPHVIZ_DEBUG && (remaining[i]->truth != NULL)) {
			fprintf(debug_graphs, "\t%u -> %u;\n", graph_label, graph_id);
			fprintf(debug_graphs, "\t%u [label=\"N=%i\\ntry=%lf\"];\n", graph_id, remaining[i]->truth->true_n, combined_error);
			part_graph_id = graph_id;
			graph_id++;
		}
		if (TREE_DEBUG) {
			printf("\tstep error: %lf, \tcombined error: %lf\n", step_error, combined_error);
		}


		if ((combined_error < best) && (step_error < MAX_SINGLE_FOM)) {
			// time to go one layer further
			scatter** new_array = build_array_no_i(remaining, remain_count, i);
			// continuing_path = NULL;

			if (GRAPHVIZ_DEBUG) {
				fprintf(debug_graphs, "%u -> ", part_graph_id);
			}

			double below = recursive_search(best, combined_error, new_energy, new_eng_uncert, loc, remaining[i], new_array, remain_count - 1, &continuing_path);

			// if (GRAPHVIZ_DEBUG) {
			// 	fprintf(debug_graphs, "} %u -> {", graph_label);
			// }

			free(new_array);
			if (below < better_find) {
				// we found a better result
				better_find = below;
				if (below < best) {
					best = below;
				}

				// begin N vs alpha work
				if (best_path != NULL) {
					fmap(best_path, free_null);
					delete_list(best_path);
					best_path = NULL;
				}
				best_path = continuing_path;
				continuing_path = NULL;
			} else {
				fmap(continuing_path, free_null);
				delete_list(continuing_path);
				continuing_path = NULL;
			}
		}
	}
	// if (TREE_DEBUG) {
	// 	printf("recursive_search: best next found N = %i\n",best_find_N);
	// }
	if (path != NULL) {
		float* best_N_return = (float*)malloc(sizeof(float) * 2);
		// this seems like it should be loc->truth
		// however, the list is misbehaving for the first call
		if (origin->truth != NULL) {
			best_N_return[0] = loc->truth->true_n;
			best_N_return[1] = loc->truth->true_eng;
		} else {
			best_N_return[0] = -1;
			best_N_return[1] = -1;
		}
		path[0] = add_to_top(best_path, best_N_return);
	}
	// if (GRAPHVIZ_DEBUG) {
	// 	fprintf(debug_graphs, "}");
	// }
	if (better_find == INFINITY) {
		better_find = best;
	}

	return better_find;
}


scatter* multi_gamma_stat_iteration(llist* history_near, llist* history_far, double sigma_per_scatter, int* best_find_array, float* alpha_eng_distro) {
	if ((history_near == NULL) || (history_far == NULL)) {
		return NULL;
	}
	// for each in history near, choose a history far and run the recursive search.
	// hang onto the best result after each point


	scatter* best_scatter = NULL;

	int len_hist_near = list_length(history_near);
	int len_hist_far = list_length(history_far);

	if (len_hist_near < 2) {
		// history1 is too short to run the recursion process.
		if (KEEP_SINGLES && len_hist_near == 1) {
			// we have a single case
			return history_near->data;
		}
		return NULL;
	}
	if (len_hist_far < 1) {
		// idk if it is possible to even get here, but we can't continue if we do
		return NULL;
	}
	// lets check the trigger: The two scatter lists must trigger at least two
	// different modules in phi
	int near_module = 0;
	llist* location = list_tail(history_near);
	while ((near_module == 0) && (location != NULL) && (location->data != NULL)) {
		near_module = phi_trigger((scatter*)(location->data), PHI_MODULES);
		location = location->up;
	}
	location = history_far;
	int far_module = 0;
	while ((far_module == 0) && (location != NULL) && (location->data != NULL)) {
		int temp = phi_trigger((scatter*)(location->data), PHI_MODULES);
		if ((abs(temp - near_module) > MODULE_SEPERATION) && (abs(temp - near_module) < (PHI_MODULES - MODULE_SEPERATION))) {
			far_module = temp;
		}
		location = location->down;
	}
	if (far_module == 0) {
		return NULL;
	}

	possible_branches += len_hist_far * factorial(len_hist_near);

	llist* hist_near_bottom = list_tail(history_near);
	// time to turn the scattering list into an array for fast access
	scatter** scatters_near = (scatter**)malloc(len_hist_near * sizeof(scatter*));
	for (int i = 0; i < len_hist_near; i++)	{
		scatters_near[i] = (scatter*)hist_near_bottom->data;
		if (TREE_DEBUG) {
			printf("mgsi: near scatter %i, N = ", i);
			if (scatters_near[i]->truth != NULL) {
				printf("%i, ", scatters_near[i]->truth->true_n);
				printf("deposit = %lf keV\n", scatters_near[i]->deposit);
			}
			}
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
	if (GENERAL_DEBUG) {
		printf("multi_gamma_stat_iteration: scatters sorted:\n");
		for (int i = 0; i < len_hist_near; i++) {
			printf("\t%lf\n", scatters_near[i]->deposit);
		}
		printf("\n");
	}

	// sorts the scatters by deposit energy, now only keep the largest LARGEST
	scatter** scatters_near_short = (scatter**)malloc((LARGEST + SKIP) * sizeof(scatter*));
	scatter** scatters_far_short = (scatter**)malloc((LARGEST + SKIP) * sizeof(scatter*));
	for (int i = 0; i < LARGEST + SKIP; i++) {
		if (i >= len_hist_near) {
			scatters_near_short[i] = NULL;
		} else {
			scatters_near_short[i] = scatters_near[i];
			if (TREE_DEBUG) {
				printf("mgsi: near scatter short %i, N = ", i);
				if (scatters_near[i]->truth != NULL) {
					printf("%i\n", scatters_near[i]->truth->true_n);
				}
			}
		}
		if (i >= len_hist_far) {
			scatters_far_short[i] = NULL;
		} else {
			scatters_far_short[i] = scatters_far[i];
		}
	}
	if (len_hist_near > LARGEST + SKIP) {
		len_hist_near = LARGEST + SKIP;
	}
	if (len_hist_far > LARGEST + SKIP) {
		len_hist_far = LARGEST + SKIP;
	}

	double best_find = sigma_per_scatter * (len_hist_near - SKIP);

	if (SCATTER_LIST_DEBUG) {
		fprintf(debug_scatter_lists, "First gamma scatters:\n");
		fprintf(debug_scatter_lists, "\tenergy order | true N |  true eng  | measured eng |      x      |      y      |      z      |\n");
		for (int i = 0; i < len_hist_near; i++) {
			if (scatters_near_short[i]->truth != NULL) {
				fprintf(debug_scatter_lists, "\t%12i | %6i | %10f | %12f |", i, scatters_near_short[i]->truth->true_n, scatters_near_short[i]->truth->true_eng, scatters_near_short[i]->deposit);
				fprintf(debug_scatter_lists, " %11f | %11f | %11f |\n", scatters_near_short[i]->loc.x, scatters_near_short[i]->loc.y, scatters_near_short[i]->loc.z);
			}
		}
		fprintf(debug_scatter_lists, "Second gamma scatters:\n");
		fprintf(debug_scatter_lists, "\tenergy order | true N |  true eng  | measured eng |      x      |      y      |      z      |\n");
		for (int i = 0; i < len_hist_far; i++) {
			if (scatters_far_short[i]->truth != NULL) {
				fprintf(debug_scatter_lists, "\t%12i | %6i | %10f | %12f |", i, scatters_far_short[i]->truth->true_n, scatters_far_short[i]->truth->true_eng, scatters_far_short[i]->deposit);
				fprintf(debug_scatter_lists, " %11f | %11f | %11f |\n", scatters_far_short[i]->loc.x, scatters_far_short[i]->loc.y, scatters_far_short[i]->loc.z);
			}
		}
	}

	llist* best_found_path = NULL;

	for (int j = 0; j < len_hist_near; j++) {

		scatter** new_array = build_array_no_i(scatters_near_short, len_hist_near, j);

		for (int i = 0; i < len_hist_far; i++) {
			llist* path = NULL;

			if (GRAPHVIZ_DEBUG) {
				fprintf(debug_graphs, "\n\ndigraph %u%i%i {\n", hist_num, i, j);
				graph_id = 0;
			}

			double try = recursive_search(best_find, 0., 511., 0., scatters_far_short[i], scatters_near_short[j], new_array, len_hist_near - 1, &path);
			if (GRAPHVIZ_DEBUG) {
				fprintf(debug_graphs, "}\n");
			}
			if (try < best_find) {
				best_scatter = scatters_near_short[j];
				best_find = try;
				if (best_found_path != NULL) {
					fmap(best_found_path, free_null);
					delete_list(best_found_path);
				}
				best_found_path = path;
			} else {
				fmap(path, free_null);
				delete_list(path);
			}
		}
		free(new_array);
	}

	// move the path into the array of N vs alpha, beta, etc.
	llist* current_location = list_head(best_found_path);
	if (best_scatter == NULL) {

	} else if (best_find_array != NULL) {
		for (int i = 0; i < FIRST_N; i++) {
			if (current_location != NULL) {
				best_find_array[i] = (int)(((float*)(current_location->data))[0]);
				alpha_eng_distro[i] = ((float*)current_location->data)[1];
				current_location = current_location->down;
			} else {
				best_find_array[i] = -1;
				alpha_eng_distro[i] = -1;
			}
		}
	}
	fmap(best_found_path, free_null);
	delete_list(best_found_path);


	free(scatters_near);
	free(scatters_far);
	free(scatters_near_short);
	free(scatters_far_short);


	return best_scatter;

}

/* double_tree_stat_iteration
 * Takes a pair of histories and sigma per scatter and finds the best solution
 * using both trees. For each i,j, the two arrays of scatters are iterated over.
 * Starting with path j->i with i as the first scatter location a FOM is calculated
 * by recursive_search. If this FOM passes the test then the second tree is searched.
 * If that FOM also passes the test then the two FOMs are combined. This combination
 * is then minimized
 */
scatter** double_tree_stat_iteration(llist* history_1, llist* history_2, double sigma_per_scatter) {
	if ((history_1 == NULL) || (history_2 == NULL)) {
		return NULL;
	}
	// for each in history near, choose a history far and run the recursive search.
	// hang onto the best result after each point


	scatter* best_scatter_1 = NULL;
	scatter* best_scatter_2 = NULL;

	int len_hist_1 = list_length(history_1);
	int len_hist_2 = list_length(history_2);

	if ((len_hist_1 < 2) || (len_hist_2 < 2)) {
		// behavior for keeping singles goes here. If we are to keep singles
		// we should return the single from any of these that have a single
		// scatter, and run the tree search on any others
		if (KEEP_SINGLES) {
			scatter** endpoints = (scatter**)calloc(2, sizeof(scatter*));
			if (len_hist_1 == 1) {
				// hist 1 has a single scatter
				endpoints[0] = history_1->data;
				if (((scatter*)(history_1->data))->truth != NULL) {
					alpha_n_1[0] = ((scatter*)(history_1->data))->truth ->true_n;
					alpha_eng_distro_n_1[0] = ((scatter*)(history_1->data))->truth ->true_eng;
				}
			}
			if (len_hist_2 == 1) {
				endpoints[1] = history_2->data;
				if (((scatter*)(history_2->data))->truth != NULL) {
					alpha_n_2[0] = ((scatter*)(history_2->data))->truth ->true_n;
					alpha_eng_distro_n_2[0] = ((scatter*)(history_2->data))->truth ->true_eng;
				}

			}
			if (endpoints[0] == NULL) {
				endpoints[0] = multi_gamma_stat_iteration(history_1, history_2, sigma_per_scatter, alpha_n_1, alpha_eng_distro_n_1);
			}
			if (endpoints[1] == NULL) {
				endpoints[1] = multi_gamma_stat_iteration(history_2, history_1, sigma_per_scatter, alpha_n_2, alpha_eng_distro_n_2);
			}
			return endpoints;
		}
		return NULL;
	}

	double best_find_1 = sigma_per_scatter * (len_hist_1 - SKIP);
	double best_find_2 = sigma_per_scatter * (len_hist_2 - SKIP);

	used_scatters_1 = (len_hist_1 - SKIP);
	used_scatters_2 = (len_hist_2 - SKIP);

	// lets check the trigger: The two scatter lists must trigger at least two
	// different modules in phi
	int near_module = 0;
	llist* location = list_tail(history_1);
	while ((near_module == 0) && (location != NULL) && (location->data != NULL)) {
		near_module = phi_trigger((scatter*)(location->data), PHI_MODULES);
		location = location->up;
	}
	location = history_2;
	int far_module = 0;
	while ((far_module == 0) && (location != NULL) && (location->data != NULL)) {
		int temp = phi_trigger((scatter*)(location->data), PHI_MODULES);
		if ((abs(temp - near_module) > MODULE_SEPERATION) && (abs(temp - near_module) < (PHI_MODULES - MODULE_SEPERATION))) {
			far_module = temp;
		}
		location = location->down;
	}
	if (far_module == 0) {
		return NULL;
	}

	possible_branches += len_hist_2 * factorial(len_hist_1);
	possible_branches += len_hist_1 * factorial(len_hist_2);

	llist* hist_1_bottom = list_tail(history_1);
	// time to turn the scattering list into an array for fast access
	scatter** scatters_1 = (scatter**)malloc(len_hist_1 * sizeof(scatter*));
	for (int i = 0; i < len_hist_1; i++)	{
		scatters_1[i] = (scatter*)hist_1_bottom->data;
		if (TREE_DEBUG) {
			printf("dtsi: near scatter %i, N = ", i);
			if (scatters_1[i]->truth != NULL) {
				printf("%i, ", scatters_1[i]->truth->true_n);
				printf("deposit = %lf keV\n", scatters_1[i]->deposit);
			}
			}
		hist_1_bottom = hist_1_bottom->up;
	}

	llist* hist_2_bottom = list_tail(history_2);
	scatter** scatters_2 = (scatter**)malloc(len_hist_2 * sizeof(scatter*));
	for (int i = 0; i < len_hist_2; i++)	{
		scatters_2[i] = (scatter*)hist_2_bottom->data;
		hist_2_bottom = hist_2_bottom->up;
	}

	scatter_quicksort(scatters_1, 0, len_hist_1 - 1);
	scatter_quicksort(scatters_2, 0, len_hist_2 - 1);
	if (GENERAL_DEBUG) {
		printf("double_tree_stat_iteration: scatters sorted:\n");
		for (int i = 0; i < len_hist_1; i++) {
			printf("\t%lf\n", scatters_1[i]->deposit);
		}
		printf("\n");
	}

	// sorts the scatters by deposit energy, now only keep the largest LARGEST
	scatter** scatters_1_short = (scatter**)malloc((LARGEST + SKIP) * sizeof(scatter*));
	scatter** scatters_2_short = (scatter**)malloc((LARGEST + SKIP) * sizeof(scatter*));
	for (int i = 0; i < LARGEST + SKIP; i++) {
		if (i >= len_hist_1) {
			scatters_1_short[i] = NULL;
		} else {
			scatters_1_short[i] = scatters_1[i];
			if (TREE_DEBUG) {
				printf("dtsi: near scatter short %i, N = ", i);
				if (scatters_1[i]->truth != NULL) {
					printf("%i\n", scatters_1[i]->truth->true_n);
				}
			}
		}
		if (i >= len_hist_2) {
			scatters_2_short[i] = NULL;
		} else {
			scatters_2_short[i] = scatters_2[i];
		}
	}
	if (len_hist_1 > LARGEST + SKIP) {
		len_hist_1 = LARGEST + SKIP;
	}
	if (len_hist_2 > LARGEST + SKIP) {
		len_hist_2 = LARGEST + SKIP;
	}

	if (SCATTER_LIST_DEBUG) {
		fprintf(debug_scatter_lists, "First gamma scatters:\n");
		fprintf(debug_scatter_lists, "\tenergy order | true N |  true eng  | measured eng |      x      |      y      |      z      |\n");
		for (int i = 0; i < len_hist_1; i++) {
			if (scatters_1_short[i]->truth != NULL) {
				fprintf(debug_scatter_lists, "\t%12i | %6i | %10f | %12f |", i, scatters_1_short[i]->truth->true_n, scatters_1_short[i]->truth->true_eng, scatters_1_short[i]->deposit);
				fprintf(debug_scatter_lists, " %11f | %11f | %11f |\n", scatters_1_short[i]->loc.x, scatters_1_short[i]->loc.y, scatters_1_short[i]->loc.z);
			}
		}
		fprintf(debug_scatter_lists, "Second gamma scatters:\n");
		fprintf(debug_scatter_lists, "\tenergy order | true N |  true eng  | measured eng |      x      |      y      |      z      |\n");
		for (int i = 0; i < len_hist_2; i++) {
			if (scatters_2_short[i]->truth != NULL) {
				fprintf(debug_scatter_lists, "\t%12i | %6i | %10f | %12f |", i, scatters_2_short[i]->truth->true_n, scatters_2_short[i]->truth->true_eng, scatters_2_short[i]->deposit);
				fprintf(debug_scatter_lists, " %11f | %11f | %11f |\n", scatters_2_short[i]->loc.x, scatters_2_short[i]->loc.y, scatters_2_short[i]->loc.z);
			}
		}
	}



	llist* best_found_path_1 = NULL;
	llist* best_found_path_2 = NULL;
	double best_combination = best_find_1 + best_find_2;

	for (int j = 0; j < len_hist_1; j++) {

		scatter** new_array_1 = build_array_no_i(scatters_1_short, len_hist_1, j);

		double reduced_second_FOM = best_find_2;

		for (int i = 0; i < len_hist_2; i++) {
			llist* path_1 = NULL;

			if (GRAPHVIZ_DEBUG) {
				fprintf(debug_graphs, "\n\ndigraph %u%i%i {\n", hist_num, i, j);
				graph_id = 0;
			}

			double try_1 = recursive_search(best_find_1, 0., 511., 0., scatters_2_short[i], scatters_1_short[j], new_array_1, len_hist_1 - 1, &path_1);
			if (GRAPHVIZ_DEBUG) {
				fprintf(debug_graphs, "}\n");
			}

			double try_2 = INFINITY;
			llist* path_2 = NULL;
			if ((try_1 < best_combination) && (try_1 < best_find_1)) {

				if (GRAPHVIZ_DEBUG) {
					fprintf(debug_graphs, "\n\ndigraph %u%i%i {\n", hist_num, j, i);
					graph_id = 0;
				}

				reduced_second_FOM = best_combination - try_1;
				if (reduced_second_FOM > best_find_2) {
					reduced_second_FOM = best_find_2;
				}

				scatter** new_array_2 = build_array_no_i(scatters_2_short, len_hist_2, i);
				try_2 = recursive_search(reduced_second_FOM, 0., 511., 0., scatters_1_short[j], scatters_2_short[i], new_array_2, len_hist_2 - 1, &path_2);
				if (GRAPHVIZ_DEBUG) {
					fprintf(debug_graphs, "}\n");
				}

				free(new_array_2);
			}


			double combination = try_1 + try_2;//(try_1 * try_1) + (try_2 * try_2);
			if ((try_1 < best_find_1) && (try_2 < reduced_second_FOM) && (combination < best_combination)) {
				best_combination = combination;
				best_scatter_1 = scatters_1_short[j];
				first_FOM = try_1;
				second_FOM = try_2;
				if (best_found_path_1 != NULL) {
					fmap(best_found_path_1, free_null);
					delete_list(best_found_path_1);
				}
				best_found_path_1 = path_1;

				best_scatter_2 = scatters_2_short[i];
				if (best_found_path_2 != NULL) {
					fmap(best_found_path_2, free_null);
					delete_list(best_found_path_2);
				}
				best_found_path_2 = path_2;
			} else {
				fmap(path_1, free_null);
				delete_list(path_1);
				fmap(path_2, free_null);
				delete_list(path_2);
			}
		}
		free(new_array_1);
	}

	// move the path into the array of N vs alpha, beta, etc.
	for (int i = 0; i < 2; i++) {
		llist* current_location = NULL;
		scatter* search_scatter = NULL;
		int *t_vs_r;
		float* E_t_vs_r;
		if (i == 0) {
			current_location = list_head(best_found_path_1);
			search_scatter = best_scatter_1;
			t_vs_r = alpha_n_1;
			E_t_vs_r = alpha_eng_distro_n_1;
		} else {
			current_location = list_head(best_found_path_2);
			search_scatter = best_scatter_2;
			t_vs_r = alpha_n_2;
			E_t_vs_r = alpha_eng_distro_n_2;
		}
		if (search_scatter != NULL) {
			for (int i = 0; i < FIRST_N; i++) {
				if (current_location != NULL) {
					t_vs_r[i] = (int)(((float*)(current_location->data))[0]);
					E_t_vs_r[i] = ((float*)current_location->data)[1];
					current_location = current_location->down;
				} else {
					t_vs_r[i] = -1;
					E_t_vs_r[i] = -1;
				}
			}
		}
	}
	fmap(best_found_path_1, free_null);
	delete_list(best_found_path_1);
	fmap(best_found_path_2, free_null);
	delete_list(best_found_path_2);


	free(scatters_1);
	free(scatters_2);
	free(scatters_1_short);
	free(scatters_2_short);

	scatter** both = (scatter**)malloc(sizeof(scatter*) * 2);
	both[0] = best_scatter_1;
	both[1] = best_scatter_2;

	return both;
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
				double photons = cur_event->energy * P_per_keV;
				add_scatter = new_scatter(scatter_loc, three_vec(0,0,0), gamma_id, cur_event->energy + deposit, cur_event->tof, sqrt(photons) / P_per_keV, spc_uncert, time_uncert_cm);
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
		double time_adjust = time_var * (time_uncert_cm / SPD_LGHT);
		// blur the time
		working_crystal->time += time_adjust;

		// calculate the new energy using the number of photons we will count
		double new_eng = (eng_var * (sqrt(working_crystal->energy * P_per_keV) / P_per_keV)) + working_crystal->energy;
		working_crystal->energy = new_eng; // set the new energy

	}

	return list_of_crystals;
}

// takes a crystal and converts it to a scatter for use by a search algorithm
scatter* crystal_to_scatter(crystal* in) {
	if (in == NULL) {
		return NULL;
	}
	scatter* new = new_scatter(in->pos, three_vec(0,0,0), in->true_gamma, in->energy, in->time, sqrt(in->energy / P_per_keV), spc_uncert, time_uncert_cm);
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

// we now have the infirsturcture for finding what scatter associated with a
// gamma is the first. Now we need to run this twice (once on each annihilation
// gamma) to get the two ends. Then from that we have the two ends and can
// run a modified version of first_scat_miss from the energy cut calculations



/* 
 * find_double_endpoints_stat
 * finds the predicted endpoints of the first two gammas in the detector history
 * using a statistical determination method. For the first step a scatter is
 * chosen from each gamma. Then recursivly a single scatter is chosen from the
 * list of remaining scatters being seached, and the probability of the deviation
 * from physical that is found is determined. If this does not exceed the current
 * best find the process repeats until the list is empty. The recursive portion
 * is completed using the recursive_search function. The energy_percent value
 * provides a minimum probability that must be cleared for the result to be
 * passed out.
 * 
 * Now uses the double_tree_stat_iteration method, moving some of the behavior 
 * done by find_endpoints_stat into double_tree_stat_iteration instead.
 */
scatter** find_double_endpoints_stat(llist* detector_history, double sigma_per_scatter) {
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
	first_id = 3;
	second_id = 2;

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
	// clear the alpha_n arrays
	for (int i = 0; i < FIRST_N; i++) {
		alpha_n_1[i] = -1;
		alpha_n_2[i] = -1;
	}
	llist* current_loc1 = list_tail(scat_list1);
	llist* current_loc2 = list_tail(scat_list2);
	for (int i = 0; i < FIRST_N; i++) {
		if ((current_loc1 != NULL) && (current_loc1->data != NULL)) {
			n_eng_distro_n_1[i] = ((scatter*)(current_loc1->data))->deposit;
			current_loc1 = current_loc1->up;
		} else {
			n_eng_distro_n_1[i] = -1;
		}
		if ((current_loc2 != NULL) && (current_loc2->data != NULL)) {
			n_eng_distro_n_2[i] = ((scatter*)(current_loc2->data))->deposit;
			current_loc2 = current_loc2->up;
		} else {
			n_eng_distro_n_2[i] = -1;
		}
	}
	if (SCATTER_LIST_DEBUG) {
		fprintf(debug_scatter_lists, "\nHistory %i:\n", ((event*)(detector_history->data))->number);
	}

	// run the actual finding of endpoint 1
	scatter** endpoints = double_tree_stat_iteration(scat_list1, scat_list2, sigma_per_scatter);

	scatter* endpoint1 = NULL;
	scatter* endpoint2 = NULL;
	if (endpoints != NULL) {
		endpoint1 = endpoints[0];
		endpoint2 = endpoints[1];
		free(endpoints);
	}
	// if we failed just give the highest energy scatter. Leads to far more bad
	// solutions but also doesn't leave anything on the table.
	if (NEVER_CUT) {
		if ((endpoint1 == NULL) && (scat_list1->data != NULL)) {
			endpoint1 = scat_list1->data;
		}
		if ((endpoint2 == NULL) && (scat_list2->data != NULL)) {
			endpoint2 = scat_list2->data;
		}
	}


	if ((endpoint1 == NULL) || (endpoint2 == NULL)) {
		fmap(scat_list1, delete_scatter);
		fmap(scat_list2, delete_scatter);
		delete_list(scat_list1);
		delete_list(scat_list2);
		missed_reconstructions += missed_reconstruction_IPS_mask;
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
	double time_delta = b->time - a->time;
	vec3d displacement = vec_scaler(ba_unit, SPD_LGHT * time_delta * 0.5);
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
		printf("history file, output file name,");
		printf(" and FOM cutoff.\n");
		// printf("");
		printf("Use a -h to get additional help information\n");
		return 1;
	}
	int binary = 0;
	if (argc > 4) {
		// go check all of the following flags!
		for (int i = 1; i < argc; i++) {
			if (!strcasecmp(argv[i], "-b")) {
				// binary flag
				binary = 1;
			}
			if (!strcasecmp(argv[i], "-E")) {
				// change energy per scatter
				if (argc > (i + 1)) {
					P_per_keV = strtod(argv[i+1], NULL);
					printf("Energy per switch set to %lf keV/switch\n", P_per_keV);
				} else {
					fprintf(stderr, "-E flag not followed by float!\n");
					return 1;
				}
			}
			if (!strcasecmp(argv[i], "-s")) {
				// change spatial resolution
				if (argc > (i + 1)) {
					spc_uncert = strtod(argv[i+1], NULL);
					printf("Spatial resolution set to %lf cm\n", spc_uncert);
				} else {
					fprintf(stderr, "-s flag not followed by float!\n");
					return 1;
				}
			}
			if (!strcasecmp(argv[i], "-t")) {
				// change time randomness
				if (argc > (i + 1)) {
					time_uncert_cm = strtod(argv[i+1], NULL);
					printf("Time resolution set to %lf cm\n", time_uncert_cm);
				} else {
					fprintf(stderr, "-t flag not followed by float!\n");
					return 1;
				}
			}
			if (!strcasecmp(argv[i], "-d")) {
				// change application of time randomness
				run_time_rand = 0;
				printf("Time randomness will not be applied\n");
			}
			if (!strcasecmp(argv[i], "-h")) {
				printf("\nHelp information for reverse kinematics:\n");
				printf("\tExpects 3 arguments, then optional flags. The three\n");
				printf("\targuments are:\n\t\t1) The DetectorTuple.pshp file containing ");
				printf("the gamma\n\t\tand electron interactions from the detector volume\n");
				printf("\t\t2) The location and name for the output files (file \n");
				printf("\t\t extensions get added automatically)\n");
				printf("\t\t3) The FOM per scatter, we have typically used 1.3\n");
				printf("\tThese arguments can then be followed by optional flags:\n");
				printf("\n\t-b -- sets code to expect tuple data in binary format\n");
				printf("\n\t-E [value in keV] -- sets keV/switch for the run\n");
				printf("\n\t-s [value in cm] -- sets spatial resolution for the run\n");
				printf("\n\t-t [value in cm] -- sets 1 sigma timing resolution IN CM\n");
				printf("\n\t-h\n\t-H -- display the help information\n");
				printf("\n\t-d -- turns off applying time randomness. This is for use\n");
				printf("\tin debugging, particularly understanding close misses\n");
				printf("\t\n");
			}
		}
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
	char* debug_file = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
	char* eng_file = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
	char* true_file;
	char* misID_file;
	char* lead_file;
	if (LOR_GROUP) {
		true_file = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
		misID_file = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
		lead_file = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
	}

	strcpy(debug_file, argv[2]);
	strcpy(lor_file, argv[2]);
	strcpy(eng_file, argv[2]);
	if (LOR_GROUP) {
		strcpy(true_file, argv[2]);
		strcpy(misID_file, argv[2]);
		strcpy(lead_file, argv[2]);
	}

	debug_file = strcat(debug_file, ".debug");
	lor_file = strcat(lor_file, ".lor");
	eng_file = strcat(eng_file, ".eng");
	if (LOR_GROUP) {
		true_file = strcat(true_file, ".true");
		misID_file = strcat(misID_file, ".misID");
	}
	

	FILE* out_in_patient = fopen(debug_file, "w");
	FILE* lor_output = fopen(lor_file, "w");
	FILE* eng_output = fopen(eng_file, "w");
	FILE* true_out;
	FILE* misID_out;
	if (LOR_GROUP) {
		true_out = fopen(true_file, "w");
		misID_out = fopen(misID_file, "w");
	}
	if ((out_in_patient == NULL) || (lor_output == NULL) || (eng_output == NULL)) {
		printf("Unable to open output file for writing\n");
		return 1;
	}

	if (GRAPHVIZ_DEBUG) {
		char* graph_file_name = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
		strcpy(graph_file_name, argv[2]);
		graph_file_name = strncat(graph_file_name, ".dot", strlen(argv[2]) + 9);
		debug_graphs = fopen(graph_file_name, "w");
	}
	if (SCATTER_LIST_DEBUG) {
		char* scatter_file_name = (char*)malloc(sizeof(char) * (strlen(argv[2]) + 10));
		strcpy(scatter_file_name, argv[2]);
		scatter_file_name = strncat(scatter_file_name, ".scatters", strlen(argv[2]) + 9);
		debug_scatter_lists = fopen(scatter_file_name, "w");
	}

	int run_in_patient = 1;
	energy_cutoff = strtod(argv[3], NULL);
	
	if ((!test_expected_energy()) || (test_vec_to_phi() != 3) || (!test_factorial()) || (test_crystal_indicies())) {
		fprintf(stderr, "tests failed, exiting\n");
		return 1;
	}

	detector_bore = shape_build(CYLINDER, bore_position, bore_dimentions, 3, -1.0);

	fprintf(out_in_patient, "history number, in patient scatter occurance, ");
	fprintf(out_in_patient, "scatters in branch 1, scatters in branch 2,");
	fprintf(out_in_patient, "R1_1, R2_1, R3_1, R4_1, R5_1, ");
	fprintf(out_in_patient, "R1_2, R2_2, R3_2, R4_2, R5_2, FOM 1, FOM 2, depth 1, depth 2\n");

	llist *in_det_hist = NULL;
	// llist *history = NULL;
	if (binary) {
		// history = load_history(in_histories, read_line_binary);
		in_det_hist = load_history(in_det_histories, read_line_binary);
	} else {
		// history = load_history(in_histories, read_line);
		in_det_hist = load_history(in_det_histories, read_line);
	}
	// some extra vars for calculating percentage of in paitent scatters
	// uint table[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	// uint table_total = 0;
	// uint numerator = 0;
	// uint denominator = 0;
	hist_num = 0;
	first_scat_hypot = 0;
	second_scat_hypot = 0;
	// uint gold_events = 0;
	// uint gold_IPS = 0;
	// uint silver_events = 0;
	// uint silver_IPS = 0;
	// uint lead_events = 0;
	// uint lead_IPS = 0;

	// begin the primary loop over all histories
	while (in_det_hist != NULL) {
		// print an update to how far we have made it
		if ((hist_num / 100000) * 100000 == hist_num) {
			printf("history number: %u\n", hist_num);
		}


		first_scat_hypot = 0;
		second_scat_hypot = 0;
		// llist* curr_loc = list_head(history);

		// does the first annhilation gamma scatter, and how much
		int wasIPS = in_patient(in_det_hist);

		// fmap(history, print_event);
		if (wasIPS && run_in_patient) {
			fprintf(out_in_patient, "%i, 1, ", ((event*)(in_det_hist->data))->number);
		} else {
			fprintf(out_in_patient, "%i, 0, ", ((event*)(in_det_hist->data))->number);
		}
		// done determining in patient scattering

		// first we need to find the location of the endpoint scatters

		// first_scat_hypot = 0; // energy differences between hypotheses
		// second_scat_hypot = 0;

		// predicted_vs_real[0] = 0;
		// predicted_vs_real[1] = 0;
		// predicted_vs_real[2] = 0;
		// predicted_vs_real[3] = 0;
		if (wasIPS) {
			missed_reconstruction_IPS_mask = 0;
		} else {
			missed_reconstruction_IPS_mask = 1;
		}
		branch_scatters_1 = 0;
		branch_scatters_2 = 0;

		scatter** endpoints = find_double_endpoints_stat(in_det_hist, energy_cutoff);
		// scatter** endpoints = find_endpoints_stat(in_det_hist, energy_cutoff);

		if (endpoints == NULL) {
			fprintf(out_in_patient, "%i, ", branch_scatters_1); // trans miss
			fprintf(out_in_patient, "%i, ", branch_scatters_2); // longi miss


		} else {
			// create the LOR
			lor* result = create_lor(endpoints[0], endpoints[1]);

			// add the in patient scatters to the number T value for scatters
			for (int i = 0; i < FIRST_N; i++) {
				if (alpha_n_1[i] > 0) {
					alpha_n_1[i] += (wasIPS & 0x1);
				}
				if (alpha_n_2[i] > 0) {
					alpha_n_2[i] += !(!(wasIPS & 0x2));
				}
			}

			if (!LOR_GROUP) {
				fprintf(lor_output, "%i, ", ((event*)(in_det_hist->data))->number);
				print_lor(lor_output, result);
				fprintf(lor_output, "\n");
			} else {
				if (CUT_IPS && (wasIPS)) {
					// we had an in patient scatter, so don't output the data
				} else {
					if ((alpha_n_1[0] == 1) && (alpha_n_2[0] == 1)) {
						fprintf(true_out, "%i, ", ((event*)(in_det_hist->data))->number);
						print_lor(true_out, result);
						fprintf(true_out, "\n");
					} else {
						fprintf(misID_out, "%i, ", ((event*)(in_det_hist->data))->number);
						print_lor(misID_out, result);
						fprintf(misID_out, "\n");
					}
				}
			}

			// now we need to find the distance by which the endpoints miss
			
			// double miss_dist_trans = first_scat_miss_transverse(result, annh_loc);
			// double miss_dist_long = first_scat_miss_longitudinal(result, annh_loc);
			// free(annh_loc);


			fprintf(out_in_patient, "%i, ", branch_scatters_1); // trans miss
			fprintf(out_in_patient, "%i, ", branch_scatters_2); // longi miss

			free_lor(result);

		}
		for (int i = 0; i < FIRST_N; i++) {
			fprintf(out_in_patient, "%i, ", alpha_n_1[i]);
			alpha_n_1[i] = -1;
		} // alpha, beta ... of scatter set 1
		for (int i = 0; i < FIRST_N; i++) {
			fprintf(out_in_patient, "%i, ", alpha_n_2[i]);
			alpha_n_2[i] = -1;
		} // alpha, beta ... of scatter set 2
		fprintf(out_in_patient, " %f, %f, ", first_FOM, second_FOM); 
		first_FOM = -1;
		second_FOM = -1;
		// print the figure of merit of the two branches
		fprintf(out_in_patient, "%i, %i\n", used_scatters_1, used_scatters_2);
		// print the depth that the search went to for both
		used_scatters_1 = -1;
		used_scatters_2 = -1;

		for (int i = 0; i < FIRST_N; i++) {
			fprintf(eng_output, "%f, ", alpha_eng_distro_n_1[i]);
			alpha_eng_distro_n_1[i] = -1;
		} // alpha, beta ... of scatter set 1
		for (int i = 0; i < FIRST_N; i++) {
			fprintf(eng_output, "%f, ", alpha_eng_distro_n_2[i]);
			alpha_eng_distro_n_2[i] = -1;
		} // alpha, beta ... of scatter set 2
		for (int i = 0; i < FIRST_N; i++) {
			fprintf(eng_output, "%f, ", n_eng_distro_n_1[i]);
			n_eng_distro_n_1[i] = -1;
		} // alpha, beta ... of scatter set 1
		for (int i = 0; i < FIRST_N; i++) {
			fprintf(eng_output, "%f, ", n_eng_distro_n_2[i]);
			n_eng_distro_n_2[i] = -1;
		} // alpha, beta ... of scatter set 2
		fprintf(eng_output, " 0\n"); // ending character

		if (endpoints != NULL) {
			delete_scatter(endpoints[0]);
			delete_scatter(endpoints[1]);
			free(endpoints);
		}

		// } else {
		// 	fprintf(out_in_patient, "%f, ", -1.); // trans miss
		// 	fprintf(out_in_patient, "%f, ", -1.); // longi miss
		// 	for (int i = 0; i < FIRST_N; i++) {
		// 		fprintf(out_in_patient, "%i, ", -1);
		// 		alpha_n_1[i] = -1;
		// 	} // alpha, beta ... of scatter set 1
		// 	for (int i = 0; i < FIRST_N; i++) {
		// 		fprintf(out_in_patient, "%i, ", -1);
		// 		alpha_n_2[i] = -1;
		// 	} // alpha, beta ... of scatter set 2
		// 	fprintf(out_in_patient, " 0\n"); // ending character
		// }


		fmap(in_det_hist, delete_event);
		delete_list(in_det_hist);
		if (binary) {
			in_det_hist = load_history(in_det_histories, read_line_binary);
		} else {
			in_det_hist = load_history(in_det_histories, read_line);
		}
		hist_num++;

	}
	// if (denominator != 0) {
	// 	printf("Percentage of in paitent scatters: ");
	// 	float percent = (float)numerator / (float)denominator;
	// 	printf("%f\n", percent);
	// }


	printf("\n\nTotal scatters: %u\nElectron path scatters%u\n",total_scatters,path_scatters);
	printf("Percentage with path: %lf\n", ((double)path_scatters/(double)total_scatters));

	printf("Number of possible branches: %llu\n", possible_branches);

	printf("Missed reconstructions (without IPS): %u\n", missed_reconstructions);
	printf("Cluster finding performace: %lf %%\n", 100. * ((float)correct_clustering / (float)total_clustering));

	fclose(out_in_patient);
	return 0;
}
