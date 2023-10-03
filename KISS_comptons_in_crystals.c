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
#define CRYSTAL_HEIGHT 1.8
#define CRYSTAL_Z_THICKNESS 0.276
char time_ordering = 0;
char cluster_energy_cut = 0;
char print_energies = 0;
FILE* energies_file = NULL;


double time_uncert_cm = 3.82; // in cm for one sigma, NOT ps or ns FWHM 300ps = 3.82
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
#define LOR_GROUP 1 // iff 1 outputs true and misID
#define CUT_IPS 0 // do not output to true/misID files if an IPS happened

//OUTDATED if 0 output all found lors, 1 only with both T=R=1, 2 
// only R!=T=1 for only 1, 3 only R!=T=1 for both 

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


typedef struct _crystal {
	double energy;
	double time;
	int index;
	uint true_n;
	int true_gamma;
	vec3d pos;
} crystal;

typedef struct _debug_data {
	double dist_1;
	double dist_2;
	char photoelectric_1;
	char photoelectric_2;
	int branch_scatters_1;
	int branch_scatters_2;
	int alpha_n_1[FIRST_N];
	int alpha_n_2[FIRST_N];
} debug;

void default_debug(debug* diagnostic) {
	if (diagnostic != NULL) {
		diagnostic->dist_1 = -1;
		diagnostic->dist_2 = -1;
		diagnostic->photoelectric_1 = 0;
		diagnostic->photoelectric_2 = 0;
		diagnostic->branch_scatters_1 = 0;
		diagnostic->branch_scatters_2 = 0;
		for (int i = 0; i < FIRST_N; i++) {
			diagnostic->alpha_n_1[i] = -1;
			diagnostic->alpha_n_2[i] = -1;
		}
	}
	return;
}
  
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
	pos.z += 200.0; // offset to avoid annoying truncation behavior
	int z_silce_crystals = crystals_in_phi * ((int)ceil((double)SCANNER_THICKNESS / crystal_r)); // crystals in a z slice
	int z_slice_index = (int)floor(pos.z / z_thickness) * z_silce_crystals; // how many z slices from 0
	// printf("pos_to: phi = %lf\n", vec_to_phi(pos));
	int phi_index = (int)floor((vec_to_phi(pos) / (PI * 2.0)) * crystals_in_phi); // how far rotated in phi
	pos.z = 0.0;
	// printf("pos_to: r = %lf\n", vec_mag(pos));
	int r_index = (int)floor((vec_mag(pos) - (double)BORE_RADIUS) / crystal_r) * crystals_in_phi;
	return phi_index + z_slice_index + r_index;
}


vec3d crystal_index_to_pos(int index, int crystals_in_phi, double crystal_r, double z_thickness) {
	int z_silce_crystals = crystals_in_phi * ((int)ceil((double)SCANNER_THICKNESS / crystal_r)); // crystals in a z slice
	double z_dist = (index / z_silce_crystals) * z_thickness + (0.5 * z_thickness); // how far in Z was this slice
	int index_no_z = index % z_silce_crystals; // mod by z crystals
	double radius = (index_no_z / crystals_in_phi) * crystal_r + BORE_RADIUS + (0.5 * crystal_r); // get the radius
	int index_no_r = index_no_z % crystals_in_phi; // mod by radius
	double phi = (((double)index_no_r + 0.5) / (double)crystals_in_phi) * (PI * 2.0); // get the phi
	vec3d cartesian;
	cartesian.z = z_dist - 200.0; // undo the z shift used to avoid a change in truncation behavior around 0
	cartesian.x = radius * cos(phi);
	cartesian.y = radius * sin(phi); // convert from cylindrical to cartesian
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
	} else {
		printf("test_crystal_indicies: passed test 1\n");
	}
	initial_pos = three_vec(40.1, 34.3, -101.1);
	// this was just a weird position to choose.
	index = pos_to_crystal_index(initial_pos, crystals_in_phi, crystal_height, z_thickness);
	final_pos = crystal_index_to_pos(index, crystals_in_phi, crystal_height, z_thickness);
	if (vec_mag(vec_sub(initial_pos, final_pos)) > 1.0) {
		// currently set to fail so that it can be checked by eye
		fprintf(stderr, "test_crystal_indicies-TEST FAILED: expected these two positions to be close:\n");
		fprintf(stderr, "\tinitial position: ");
		vec_print(initial_pos, stderr);
		fprintf(stderr, "\n\tend position: ");
		vec_print(final_pos, stderr);
		fprintf(stderr, "\n\tindex %i\n", index);
		pass += 1;
	} else {
		printf("test_crystal_indicies: passed test 2\n");
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

// selects a scatter randomly
scatter* random_selection(scatter** cluster, int length, double arg) {
    if (cluster == NULL) {
        return NULL;
    }
    int index = rand() % length;
    return cluster[index];
}

// selects scatter by the the highest energy
scatter* energy_max_selection(scatter** cluster, int length, double arg) {
    if (cluster == NULL) {
        return NULL;
    }
    return cluster[0]; // cluster should be sorted by energy
}

// selects scatter by the the highest energy
scatter* energy_min_selection(scatter** cluster, int length, double arg) {
    if (cluster == NULL) {
        return NULL;
    }
    return cluster[length - 1]; // cluster should be sorted by energy
}

scatter* time_selection(scatter** cluster, int length, double arg) {
    if (cluster == NULL) {
        return NULL;
    }
    int best_index = 0;
    double best_time = cluster[0]->time;
    for (int i = 1; i < length; i++) {
        if (cluster[i]->time < best_time) {
            best_index = i;
            best_time = cluster[i]->time;
        }
    }
    return cluster[best_index];
}

scatter* radius_selection(scatter** cluster, int length, double arg) {
    if (cluster == NULL) {
        return NULL;
    }
    int best_index = 0;
    double closest = INFINITY;
    for (int i = 0; i < length; i++) {
        vec3d pos = cluster[i]->loc;
        pos.z = 0;
        double distance = vec_mag(pos);
        if (distance < closest) {
            best_index = i;
            closest = distance;
        }
    }
    return cluster[best_index];
}

// returns any scatter with over the energy cut in energy
scatter* pe_selection(scatter** cluster, int length, double energy_cut) {
    if (cluster == NULL) {
        return NULL;
    }
	for (int i = 0; i < length; i++) {
		if (fabs(cluster[i]->deposit - ELECTRON_MASS) < energy_cut * ELECTRON_MASS) {
			return cluster[i];
		}
	}
	return NULL;
}

/* single_cluster_call
 * Takes a single cluster with more than one scatter and selects a single scatter as the "first". 
 * It does not check that it has more than one scatter so that must be enforced outside of this
 * function. The selection method is passed to this function. If cluster energy cuts are to be
 * run then this function can return NULL in event the cluster fails the cut.
 */
scatter* single_cluster_call(llist* cluster, double energy_cut, debug* diagnostic, int index_cluster, scatter* (*selection_method)(scatter**, int, double)) {
    if (cluster == NULL) {
        return NULL;
    }
	llist* cluster_bottom = list_tail(cluster);
    int len_cluster = list_length(cluster);
	// time to turn the scattering list into an array for fast access
	scatter** scatters = (scatter**)malloc(len_cluster * sizeof(scatter*));
	for (int i = 0; i < len_cluster; i++)	{
		scatters[i] = (scatter*)cluster_bottom->data;
		cluster_bottom = cluster_bottom->up;
	}
    if (diagnostic != NULL) {
        if (index_cluster == 0) {
            diagnostic->dist_1 = vec_dist(scatters[0]->loc, scatters[1]->loc);
        } else {
            diagnostic->dist_2 = vec_dist(scatters[0]->loc, scatters[1]->loc);
        }
    }
	scatter_quicksort(scatters, 0, len_cluster - 1);

    if (SCATTER_LIST_DEBUG) {
        if (index_cluster == 0) {
            fprintf(debug_scatter_lists, "First ");
        } else {
            fprintf(debug_scatter_lists, "Second ");
        }
		fprintf(debug_scatter_lists, "gamma scatters:\n");
		fprintf(debug_scatter_lists, "\tenergy order | true N |  true eng  | measured eng |      x      |      y      |      z      |\n");
		for (int i = 0; i < len_cluster; i++) {
			if (scatters[i]->truth != NULL) {
				fprintf(debug_scatter_lists, "\t%12i | %6i | %10f | %12f |", i, scatters[i]->truth->true_n, scatters[i]->truth->true_eng, scatters[i]->deposit);
				fprintf(debug_scatter_lists, " %11f | %11f | %11f |\n", scatters[i]->loc.x, scatters[i]->loc.y, scatters[i]->loc.z);
			}
		}
	}

    // if cluster energy cut is set check to see if we should reject this instance.
    if (cluster_energy_cut) {
        double total_energy = 0;
        for (int i = 0; i < len_cluster; i++) {
            total_energy += scatters[i]->deposit;
        }
        if (fabs(total_energy - ELECTRON_MASS) > energy_cut * ELECTRON_MASS) {
            free(scatters);
            return NULL;
        }
    }

    scatter* chosen = selection_method(scatters, len_cluster, energy_cut);
    free(scatters);
    if ((diagnostic != NULL) && (chosen != NULL) && (chosen->truth != NULL)) {
        if (index_cluster == 0) {
            diagnostic->alpha_n_1[0] = chosen->truth->true_n;
        } else {
            diagnostic->alpha_n_2[0] = chosen->truth->true_n;
        }
    }
    return chosen;

}


/* double_histories_search
 * Takes a pair of clusters and an energy cut and finds the a single LOR solution
 * using both trees. If one or more history has only a single scatter then it checks 
 * that against the energy cut to see if it can be accepted as a photoelectric interaction.
 * For clusters with more than one scatter the given method is used for choosing a scatter.
 */
scatter** double_histories_search(llist* cluster_1, llist* cluster_2, double energy_cut, scatter* (*selection_method)(scatter**, int, double), debug* diagnostic) {
	if ((cluster_1 == NULL) || (cluster_2 == NULL)) {
		return NULL;
	}
	diagnostic->photoelectric_1 = 0;
	diagnostic->photoelectric_2 = 0;

	// for each in history near, choose a history far and run the recursive search.
	// hang onto the best result after each point

    scatter** endpoints = (scatter**)calloc(2, sizeof(scatter*));

	int len_cluster_1 = list_length(cluster_1);
	int len_cluster_2 = list_length(cluster_2);

	if ((len_cluster_1 < 2) || (len_cluster_2 < 2)) {
		// behavior for keeping singles goes here. If we are to keep singles
		// we should return the single from any of these that have a single
		// scatter, and run the tree search on any others
		if (len_cluster_1 == 1) {
			// hist 1 has a single scatter
			if (fabs(((scatter*)(cluster_1->data))->deposit - ELECTRON_MASS) < (energy_cut * ELECTRON_MASS)) {
				endpoints[0] = ((scatter*)(cluster_1->data));
				diagnostic->alpha_n_1[0] = ((scatter*)(cluster_1->data))->truth ->true_n;
				diagnostic->photoelectric_1 = 1;
			}
		}
		if (len_cluster_2 == 1) {
			if (fabs(((scatter*)(cluster_2->data))->deposit - ELECTRON_MASS) < (energy_cut * ELECTRON_MASS)) {
				endpoints[1] = ((scatter*)(cluster_2->data));
				diagnostic->alpha_n_2[0] = ((scatter*)(cluster_2->data))->truth ->true_n;
				diagnostic->photoelectric_2 = 1;
			}
		}
        if ((endpoints[0] != NULL) && (endpoints[1] != NULL)) {
		    return endpoints;
        }
	}

    if ((endpoints[0] == NULL) && (len_cluster_1 > 1)) {
        endpoints[0] = single_cluster_call(cluster_1, energy_cut, diagnostic, 0, selection_method);
    }
    if ((endpoints[1] == NULL) && (len_cluster_2 > 1)) {
        endpoints[1] = single_cluster_call(cluster_2, energy_cut, diagnostic, 1, selection_method);
    }
    if ((endpoints[0] == NULL) || (endpoints[1] == NULL)) {
        free(endpoints);
        return NULL;
    }

	return endpoints;
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
	llist* printing_list = list_head(list_of_crystals);
	while (printing_list != NULL) {
		fprintf(energies_file, "%lf\n", ((crystal*)(printing_list->data))->energy);
		printing_list = printing_list->down;
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
scatter** find_double_endpoints_stat(llist* detector_history, double energy_cut, scatter* (*selection_method)(scatter**, int, double), debug* diagnostics) {
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

	if (SCATTER_LIST_DEBUG) {
		fprintf(debug_scatter_lists, "\nHistory %i:\n", ((event*)(detector_history->data))->number);
	}

	diagnostics->branch_scatters_1 = list_length(scat_list1);
	diagnostics->branch_scatters_2 = list_length(scat_list2);
	// run the actual finding of endpoint 1
	scatter** endpoints = double_histories_search(scat_list1, scat_list2, energy_cut, selection_method, diagnostics);

	scatter* endpoint1 = NULL;
	scatter* endpoint2 = NULL;
	if (endpoints != NULL) {
		endpoint1 = endpoints[0];
		endpoint2 = endpoints[1];
		free(endpoints);
	}
	// if we failed just give the highest energy scatter. Leads to far more bad
	// solutions but also doesn't leave anything on the table.
	// if (NEVER_CUT) {
	// 	if ((endpoint1 == NULL) && (scat_list1->data != NULL)) {
	// 		endpoint1 = scat_list1->data;
	// 	}
	// 	if ((endpoint2 == NULL) && (scat_list2->data != NULL)) {
	// 		endpoint2 = scat_list2->data;
	// 	}
	// }


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

	new->transverse_uncert = CRYSTAL_Z_THICKNESS * 0.5;
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
	int exit = 0;
	if (argc < 5) {
		printf("Unable to run. Expected at least 4 arguments got %i.\n", argc - 1);
		printf("Expects an input detector volume ");
		printf("history file, output file name,");
		printf(" FOM cutoff, and energy cut\n");
		// printf("");
		printf("Use a -h to get additional help information\n");
		exit = 1;
	}
	int binary = 0;
    scatter* (*selection_method)(scatter**, int, double)  = random_selection;
	if (argc > 1) {
		// go check all of the following flags!
		for (int i = 1; i < argc; i++) {
			if (!strcasecmp(argv[i], "-b")) {
				// binary flag
				binary = 1;
			}
			if (!strcasecmp(argv[i], "-Ex")) {
				// set to selection by energy
                selection_method = energy_max_selection;
			}
			if (!strcasecmp(argv[i], "-Ei")) {
				// set to selection by energy
                selection_method = energy_min_selection;
			}
			if (!strcasecmp(argv[i], "-r")) {
				// set to selection by radius
				selection_method = radius_selection;
			}
			if (!strcasecmp(argv[i], "-t")) {
				// set to selection by time
				selection_method = time_selection;
			}
			if (!strcasecmp(argv[i], "-p")) {
				// set to selection by photoelectric
				selection_method = pe_selection;
			}
			if (!strcasecmp(argv[i], "-a")) {
				// tell program to print out crystal energies
				printf("printing crystal energies\n");
				print_energies = 1;
			}
			if (!strcasecmp(argv[i], "-d")) {
				// change application of time randomness
				run_time_rand = 0;
				printf("Time randomness will not be applied\n");
			}
			if (!strcasecmp(argv[i], "-c")) {
				// set cluster energy cut application
				cluster_energy_cut = 1;
                printf("Cluster energy cut will be applied\n");
			}
			if (!strcasecmp(argv[i], "-h")) {
				printf("\nHelp information for reverse kinematics:\n");
				printf("\tExpects 4 arguments, then optional flags. The four\n");
				printf("\targuments are:\n\t\t1) The DetectorTuple.pshp file containing ");
				printf("the gamma\n\t\tand electron interactions from the detector volume\n");
				printf("\t\t2) The location and name for the output files (file \n");
				printf("\t\t extensions get added automatically)\n");
				printf("\t\t3) The FOM per scatter, we have typically used 1.3\n");
				printf("\t\t4) The energy cut (typically 0.1, giving +-10%%)\n");
				printf("\tThese arguments can then be followed by optional flags:\n");
				printf("\n\t-b  -- sets code to expect tuple data in binary format\n");
				printf("\n\t-r  -- sets to radius based Compton selection\n");
				printf("\n\t-Ex -- sets to maximum energy based Compton selection\n");
				printf("\n\t-Ei -- sets to minimum energy based Compton selection\n");
				printf("\n\t-t  -- sets to time based Compton selection\n");
				printf("\n\t-p  -- sets to photoelectric based selection\n");
				printf("\n\t-c  -- requires clusters to pass an energy cut\n");
				printf("\n\t-a  -- prints the energies of every crystal to an output file\n");
				printf("\n\t-h or -H -- display the help information\n");
				printf("\n\t-d  -- turns off applying time randomness. This is for use\n");
				printf("\tin debugging, particularly understanding close misses\n");
				printf("\t\n");
				if (exit) {
					return 1;
				}
			}
		}
	}
    printf("Compton selection method set to ");
    if (selection_method == random_selection) {
        printf("random\n");
    } else if (selection_method == energy_max_selection) {
        printf("max energy\n");
    } else if (selection_method == energy_min_selection) {
        printf("min energy\n");
    } else if (selection_method == radius_selection) {
        printf("radius\n");
    } else if (selection_method == pe_selection) {
        printf("photoelectric\n");
    } else if (selection_method == time_selection) {
        printf("time\n");
    } else {
        fprintf(stderr, "selection method does not match valid methods!\n");
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
	energies_file = fopen(eng_file, "w");
	FILE* true_out;
	FILE* misID_out;
	if (LOR_GROUP) {
		true_out = fopen(true_file, "w");
		misID_out = fopen(misID_file, "w");
	}
	if ((out_in_patient == NULL) || (lor_output == NULL)) {// || (eng_output == NULL)) {
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
	double energy_cut = strtod(argv[3], NULL);
	
	if ((test_vec_to_phi() != 3) || (!test_factorial()) || (test_crystal_indicies())) {
		fprintf(stderr, "tests failed, exiting\n");
		return 1;
	}

	fprintf(out_in_patient, "history number, in patient scatter occurance, ");
	fprintf(out_in_patient, "scatters in branch 1, scatters in branch 2,");
	fprintf(out_in_patient, "R1_1, R2_1, R3_1, R4_1, R5_1, ");
	fprintf(out_in_patient, "R1_2, R2_2, R3_2, R4_2, R5_2, FOM 1, FOM 2, depth 1, depth 2,");
	fprintf(out_in_patient, " photoelectric 1?, photoelectric 2?, first two spacing 1, first two spacing 2\n");

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
		debug diagnostic;
		default_debug(&diagnostic);

		scatter** endpoints = find_double_endpoints_stat(in_det_hist, energy_cut, selection_method, &diagnostic);
		// scatter** endpoints = find_endpoints_stat(in_det_hist, energy_cutoff);

		if (endpoints == NULL) {
			fprintf(out_in_patient, "%i, ", diagnostic.branch_scatters_1); // scatters in branch 1
			fprintf(out_in_patient, "%i, ", diagnostic.branch_scatters_2); // scatters in branch 2


		} else {
			// create the LOR
			lor* result = create_lor(endpoints[0], endpoints[1]);

			// add the in patient scatters to the number T value for scatters
			for (int i = 0; i < FIRST_N; i++) {
				if (diagnostic.alpha_n_1[i] > 0) {
					diagnostic.alpha_n_1[i] += (wasIPS & 0x1);
				}
				if (diagnostic.alpha_n_2[i] > 0) {
					diagnostic.alpha_n_2[i] += !(!(wasIPS & 0x2));
				}
			}

			fprintf(lor_output, "%i, ", ((event*)(in_det_hist->data))->number);
			print_lor(lor_output, result);
			fprintf(lor_output, "\n");
			
			if (LOR_GROUP) {
				if (CUT_IPS && (wasIPS)) {
					// we had an in patient scatter, so don't output the data
				} else {
					if ((diagnostic.alpha_n_1[0] == 1) && (diagnostic.alpha_n_2[0] == 1)) {
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
			fprintf(out_in_patient, "%i, ", diagnostic.alpha_n_1[i]);
		} // alpha, beta ... of scatter set 1
		for (int i = 0; i < FIRST_N; i++) {
			fprintf(out_in_patient, "%i, ", diagnostic.alpha_n_2[i]);
		} // alpha, beta ... of scatter set 2
		fprintf(out_in_patient, " %f, %f, ", first_FOM, second_FOM); 
		first_FOM = -1;
		second_FOM = -1;
		// print the figure of merit of the two branches
		fprintf(out_in_patient, "%i, %i, ", used_scatters_1, used_scatters_2);
		// print the depth that the search went to for both
		used_scatters_1 = -1;
		used_scatters_2 = -1;
		fprintf(out_in_patient, "%i, %i,", diagnostic.photoelectric_1, diagnostic.photoelectric_2);
		fprintf(out_in_patient, " %f, %f\n", diagnostic.dist_1, diagnostic.dist_2);




		if (endpoints != NULL) {
			delete_scatter(endpoints[0]);
			delete_scatter(endpoints[1]);
			free(endpoints);
		}



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
