#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"
#include <pthread.h>
#include "ray_trace.h"
#include "llist.h"

#define MAX_THREADS 4

#define NORMAL_TABLE_SIZE 10000

pthread_t tid[MAX_THREADS];
pthread_mutex_t volume_lock;

struct _add_lor_union {
	render* universe;
	lor* lor;
	pthread_mutex_t lock;
	pthread_cond_t cond;
	char flag;
	char die;
};


struct _add_lor_union working_lors[MAX_THREADS];
pthread_cond_t main_hold;
pthread_mutex_t main_lock;

geometry* objects = NULL;

double (*long_func)(double) = NULL;

double long_adjust = 1.0;
double trans_adjust = 1.0;

int long_binary = 0;
int trans_binary = 0;

unsigned long long int gauss_lookup_success = 0;
unsigned long long int gauss_lookup_fail = 0;
long int not_added = 0;
long int iterated = 0;

typedef struct _intvec {
	int a;
	int b;
	int c;
} int_vec;


// reads a line from a lor file. BE WARNED: it does not care about line endings,
// if your file does not always have the same number of columns it will break.
// expects a lor file with: hist_num, c_x, c_y, c_z, d_x, d_y, d_z, sigma_l, sigma_t
lor* read_lor(FILE* input) {
	if (input == NULL) {
		return NULL;
	}

	double center_x;
	double center_y;
	double center_z;
	double dir_x;
	double dir_y;
	double dir_z;
	double longitudinal;
	double transverse;

	int worked;

	worked = fscanf(input, "%*i,%lf,", &center_x);
	worked = fscanf(input, "%lf,", &center_y);
	worked = fscanf(input, "%lf,", &center_z);
	worked = fscanf(input, "%lf,", &dir_x);
	worked = fscanf(input, "%lf,", &dir_y);
	worked = fscanf(input, "%lf,", &dir_z);
	worked = fscanf(input, "%lf,", &longitudinal);
	worked = fscanf(input, "%lf", &transverse);

	if (worked == EOF) {
		return NULL;
	}
	
	// make a new event to be passed out
	lor* new = (lor*)malloc(sizeof(lor));
	if (new == NULL) {
		return NULL;
	}
	
	new->center = three_vec(center_x, center_y, center_z);
	vec3d* cleanup = three_vec(dir_x, dir_y, dir_z);
	new->dir = vec_norm(cleanup);
	free(cleanup);
	new->long_uncert = longitudinal * long_adjust;
	// fprintf(stdout, "%lf\n", new->long_uncert);
	new->transverse_uncert = transverse * trans_adjust;


	return new;
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

void print_lor(FILE* output, lor* lor) {
	fprintf(output, "%f, %f, %f,", lor->center->x, lor->center->y, lor->center->z);
	fprintf(output,  " %f, %f, %f,", lor->dir->x, lor->dir->y, lor->dir->z);
	fprintf(output, " %f, %f", lor->long_uncert, lor->transverse_uncert);
}

// takes a render (for the array), a double to add to it, and the location to
// do so. The indexes array should be of length 3, an x, y, and z coord.
// does basic addition of doubles
void add_double(render* u, double x, int* indexes, double atten) {
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += x * atten;
	return;
}

// takes a render (for the array), a double to add to it, and the location to
// do so. The indexes array should be of length 3, an x, y, and z coord.
// does addition of doubles and then also a multiplicative combination with
// log reduction
void add_mult_log(render* u, double x, int* indexes, double arg) {
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += x;
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += log((u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]]) * x + 1.0);
	return;
}

// takes a render (for the array), a double to add to it, and the location to
// do so. The indexes array should be of length 3, an x, y, and z coord.
// does addition of doubles and then also a multiplicative combination with
// log reduction
void add_mult(render* u, double x, int* indexes, double arg) {
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += x;
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += ((u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]]) * (x + 1));
	return;
}

// takes a render (for the arry), a double to add to it, and the location to do
// so. The indexes should be of length 3, and x, y, and z coord. Adds a 1 to the
// volume no matter what. You read that right. Its binary addition. If you call
// it, it will add.
void add_binary(render* u, double x, int* indexes, double atten) {
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += 1.0 * atten;
	return;
}

// takes a render (for the arry), a double to add to it, and the location to do
// so. The indexes should be of length 3, and x, y, and z coord. Adds the log of
// the value given
void add_log(render* u, double x, int* indexes, double arg) {
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += log(x + 1);
	return;}

// takes two corner defintions and rebuilds them so that they have one of all
// low values and one of all high values. Makes iteration simpler.
// All operations done in situe with given locations
void low_high_corner(vec3d* corner1, vec3d* corner2, vec3d* low_corner, vec3d* high_corner) {
	if ((corner1 == NULL) || (corner2 == NULL) || (low_corner == NULL) || (high_corner == NULL)) {
		return;
	}
	if (corner1->x < corner2->x) {
		low_corner->x = corner1->x;
		high_corner->x = corner2->x;
	} else {
		low_corner->x = corner2->x;
		high_corner->x = corner1->x;
	}
	if (corner1->y < corner2->y) {
		low_corner->y = corner1->y;
		high_corner->y = corner2->y;
	} else {
		low_corner->y = corner2->y;
		high_corner->y = corner1->y;
	}	
	if (corner1->z < corner2->z) {
		low_corner->z = corner1->z;
		high_corner->z = corner2->z;
	} else {
		low_corner->z = corner2->z;
		high_corner->z = corner1->z;
	}
}

render* read_render_def(FILE* input) {
	if (input == NULL) {
		fprintf(stderr, "read_render_def: null pointer given to render definitions");
		return NULL;
	}
	double x, y, z;
	int worked;

	worked = fscanf(input, "%lf,%lf,%lf,", &x, &y, &z); // render cube corner 1
	vec3d* corner1 = three_vec(x, y, z);
	worked = fscanf(input, "%lf,%lf,%lf", &x, &y, &z); // render cube corner 2
	vec3d* corner2 = three_vec(x, y, z);
	vec3d* low_corner = three_vec(0,0,0);
	vec3d* high_corner = three_vec(0,0,0);
	low_high_corner(corner1, corner2, low_corner, high_corner);
	free(corner1);
	free(corner2);
	int a, b, c;
	worked = fscanf(input, "%i,%i,%i", &a, &b, &c); // gets the number of voxels in each dimension
	
	double cutoff;
	worked = fscanf(input, "%lf", &cutoff);

	char method[50];
	worked = fscanf(input, "%s", method);

	if (worked == EOF) {
		return NULL;
	}

	// define the conversion from length space to voxel space
	double x_convert = a / (high_corner->x - low_corner->x);
	double y_convert = b / (high_corner->y - low_corner->y);
	double z_convert = c / (high_corner->z - low_corner->z);


	render* new = (render*)malloc(sizeof(render));
	double* volume = (double*)malloc(sizeof(double) * a * b * c);
	// allocate the volume for all of our rendering

	new->volume = volume;
	new->dimensions[0] = a;
	new->dimensions[1] = b;
	new->dimensions[2] = c;
	new->conversion = three_vec(x_convert, y_convert, z_convert);
	new->least_corner = low_corner;
	new->max_corner = high_corner;
	new->cutoff = cutoff;

	// printf("%s\n", method);

	if (strncasecmp("addition", method, 8) == 0) {
		// fprintf(stderr, "combine using addition\n");
		new->combiner = add_double;
	} else if (strncasecmp("add_multiplication_log", method, 22) == 0) {
		// fprintf(stderr, "combine using addition + log(multiplication)\n");
		new->combiner = add_mult_log;
	} else if (strncasecmp("add_multiplication", method, 18) == 0) {
		new->combiner = add_mult;
	} else if (strncasecmp("add_binary", method, 10) == 0) {
		new->combiner = add_binary;
	} else if (strncasecmp("add_log", method, 7) == 0) {
		new->combiner = add_log;
	} else {
		// for if no coorect value was given (or I messed up checking them)
		fprintf(stderr, "WARN: Mo addition method matched, using addition.\n");
		new->combiner = add_double;
	}



	// add stuff to handle various choices of conversion here when needed.

	return new;
}

// reads a file and gets all of the geometry definitions from it
// a shape is defined as: type, center_x, center_y, center_z, dim_1, dim_2,
// dim_3, axis, attenuation
geometry* read_geometry(FILE* input) {
	if (input == NULL) {
		fprintf(stderr, "WARN: no geometry found, continuing without attenuation\n");
		return NULL;
	}
	char type[64];
	int worked;
	llist* llist_geo = NULL;
	int cont = 1;
	int geo_size = 0;

	while (cont) {
		worked = fscanf(input, "%s,", type);
		float* posit = (float*)malloc(sizeof(float) * 3);
		worked = fscanf(input, "%f, %f, %f,", &posit[0], &posit[1], &posit[2]);
		float* dims = (float*)malloc(sizeof(float) * 3);
		worked = fscanf(input, "%f, %f, %f,", &dims[0], &dims[1], &dims[2]);
		int axis;
		float attenuation;
		worked = fscanf(input, "%i, %f", &axis, &attenuation);
		int type_int;
		if (strncasecmp(type, "cylinder", 8) == 0) {
			type_int = CYLINDER;
		} else if (strncasecmp(type, "prism", 5) == 0) {
			type_int = REC_PRISM;
		} else if (strncasecmp(type, "sphere", 6) == 0) {
			type_int = SPHERE;
		} else {
			type_int = 0;
			fprintf(stderr, "ERROR: shape has non-valid name\n");
		}
		if (worked != EOF) {
			shape* new = shape_build(type_int, posit, dims, axis, attenuation);
			llist_geo = add_to_bottom(llist_geo, new);
			geo_size++;
		} else {
			cont = 0;
		}
	}
	// now that we have gone through all of the shapes to be added, make them
	// into a single array
	llist* working_list = list_head(llist_geo);
	shape** array = (shape**)malloc(sizeof(shape*) * geo_size);
	for (int i = 0; (i < geo_size) && (working_list != NULL); i++) {
		array[i] = (shape*)(working_list->data);
		working_list->data = NULL;
		working_list = working_list->down;
	}
	delete_list(llist_geo);
	geometry* new = (geometry*)malloc(sizeof(geometry));
	new->size = geo_size;
	new->geo = array;
	return new;
}

/* 
 * space_to_vox
 * converts from spacial coordiantes (typically in cm) to the voxel coordinates.
 * This conversion does not do the final change to integeriaztion. The values
 * returns will likely still need to be rounded or similar to define the final
 * values.
 */
vec3d* space_to_vox(vec3d* a, render* universe) {
	if (a == NULL || universe == NULL) {
		return NULL;
	}
	// first handle the offset
	vec3d* offsetless = vec_sub(a, universe->least_corner);
	// now multiply by conversion
	offsetless->x = offsetless->x * universe->conversion->x;
	offsetless->y = offsetless->y * universe->conversion->y;
	offsetless->z = offsetless->z * universe->conversion->z;

	return offsetless;
}

/* 
 * vox_to_space
 * converts from voxel coordinates (typically orginating as integer indices) to
 * the spacial coordinates the LORs are given in (typically in cm)
 */
vec3d* vox_to_space(vec3d* a, render* universe) {
	if (a == NULL || universe == NULL) {
		return NULL;
	}
	double x = a->x / universe->conversion->x;
	double y = a->y / universe->conversion->y;
	double z = a->z / universe->conversion->z;

	vec3d* multiplied = three_vec(x, y, z);
	vec3d* offset = vec_add(multiplied, universe->least_corner);
	free(multiplied);

	return offset;
}

// ceilings a vector
void vec_ceil(vec3d *a) {
	a->x = ceil(a->x);
	a->y = ceil(a->y);
	a->z = ceil(a->z);
}

// floors a vector
void vec_floor(vec3d *a) {
	a->x = floor(a->x);
	a->y = floor(a->y);
	a->z = floor(a->z);
}

/*
 * coord_transfer
 * takes an array of 3 doubles and an array of 3 integers. It swaps the doubles
 * based on the numbers in the integers. The values in integers must range from
 * 0 to 2 and not repeat.
 * ex:
 * 	transfer[a, b, c], spots[i,j,k];
 * 		i=2,j=0,k=1;
 * 	rearrange transfer to [b,c,a]
 * The rearrangement is done in place. The passed array is modified
 */
void coord_transfer(double* transfer, int* spots) {
	if ((spots[0] == spots[1]) || (spots[0] == spots[2]) || (spots[1] == spots[2])) {
		fprintf(stderr, "coord_transfer: non-unique transfer");
		return;
	}
	for (int i = 0; i < 3; i++) {
		if ((spots[i] >= 0) && (spots[i] < 3)) {
			fprintf(stderr, "coord_transfer: transfer to location out of array");
			return;
		}
	}
	double x = transfer[0];
	double y = transfer[1];
	double z = transfer[2];
	transfer[spots[0]] = x;
	transfer[spots[1]] = y;
	transfer[spots[2]] = z;
}

/*
 * pt_to_lor_long
 * Gives the longitudinal distance from center of the LOR to the given point.
 * Done using a simple dot product using the fact the direction vector is
 * a normalized vector
 */
double pt_to_lor_long(lor* line, vec3d* pt) {
	if (line == NULL || pt == NULL) {
		return -1;
	}
	vec3d* to_point = vec_sub(pt, line->center);
	double dist = vec_dot(to_point, line->dir);
	free(to_point);
	return fabs(dist);
}

/*
 * pt_to_lor_trans
 * Gives the transverse distance from the center of the lor to the given point.
 * Done using a cross product and the fact the direction vector is normalized
 */
double pt_to_lor_trans(lor* line, vec3d* pt) {
	if (line == NULL || pt == NULL) {
		return -1;
	}
	vec3d* to_point = vec_sub(pt, line->center);
	vec3d* cross = vec_cross(to_point, line->dir);
	double dist = vec_mag(cross);
	free(to_point);
	free(cross);
	return fabs(dist);
}

/*
 * lor_box_offset
 * Takes a lor and make the offset for half of the processing box. The point
 * is defined as being at cutoff * sigma_l + (cutoff * sigma_t) in the x, y, and
 * z direction that the lor points in). This guarrentees that the entire volume
 * of the lor within the cutoff band is contained in the volume.
 */
vec3d* lor_box_offset(lor* lor, double cutoff_sigma) {
	vec3d* longitudinal = vec_scaler(lor->dir, lor->long_uncert * cutoff_sigma);
	double x = (lor->dir->x > 0) - (lor->dir->x < 0); // 1 if heading in x dir, -1 if in -x dir 0 if 0
	double y = (lor->dir->y > 0) - (lor->dir->y < 0);
	double z = (lor->dir->z > 0) - (lor->dir->z < 0);
	vec3d* trans_dir = three_vec(x, y, z);
	vec3d* trans_vec = vec_scaler(trans_dir, cutoff_sigma * lor->transverse_uncert);
	vec3d* total = vec_add(longitudinal, trans_vec);
	free(longitudinal);
	free(trans_dir);
	free(trans_vec);
	return total;
}


// returns the value of the normal distribution. Normalized to be 1 at x=0.
// x is value/sigma in a gaussian distribution. i.e. it already has been divided
// by sigma
double centered_normal(double x) {
	// static double table[NORMAL_TABLE_SIZE];
	// static int first_run = 0;
	// static double conversion;
	// // if this is our first run of the centered normal function we need to build
	// // the truth table. We will build it from 0 to the cutoff of the renderer
	// if (!first_run) {
	// 	conversion = NORMAL_TABLE_SIZE / master_copy->cutoff;
	// 	for (int i = 0; i < NORMAL_TABLE_SIZE; i++) {
	// 		double fraction = (double)i / conversion;
	// 		table[i] = exp(-0.5 * (fraction * fraction));
	// 	}
	// 	first_run = 1;
	// }
	double fraction = x;
	// double fraction = fabs(x);
	// if (fraction >= master_copy->cutoff) {
		// outside the volume of table
		// gauss_lookup_fail++;
	double exponent = -0.5 * (fraction * fraction);
	return exp(exponent);
	// }
	// just look it up!
	// gauss_lookup_success++;
	// return table[(int)(fraction * conversion)];
}
double centered_binary(double x) {
	return 1.0;
}

double centered_laplace(double x) {
	double fraction = x;
	double exponent = fabs(fraction);
	return exp(exponent);
}

/* 
 * takes a lor and a geometry and finds out what the probability of attenuation
 * was for the given lor. To do this it projects a ray forward along the lor
 * with propagate and one backwards along the lor. The combined distances times
 * attenuation coefficents are then added together to get the normalization
 */
double atten_correction(lor* lor) {
	if ((lor == NULL) || (objects == NULL)) {
		return 1.0;
	}
	vec3d* dir_back = vec_sub(NULL, lor->dir);
	ray* ray_1 = ray_build(vec_copy(lor->center), vec_norm(lor->dir));
	ray* ray_2 = ray_build(vec_copy(lor->center), vec_norm(dir_back));
	free(dir_back);
	double atten = propagate(ray_1, objects);
	atten += propagate(ray_2, objects);
	ray_free(ray_1);
	ray_free(ray_2);
	return exp(atten);
}

/*
 * add_lor
 * takes a universe and line of response. Adds the line of response to the
 * universe. The area it iterates over is the cube that fits the entire line of
 * response for the given cutoff within it. The iteration volume is a single
 * cube. See add_lor_plane for the version that makes a new iteration plane for
 * each iteration of the third axis
 */
void add_lor(render* universe, lor* lor) {
	if (universe == NULL || lor == NULL) {
		return;
	}
	// get the attenuation correction value
	double attenuation = atten_correction(lor);

	// print_lor(stdout, lor);
	// printf("\n");
	// first define the volume in which we will be operating
	vec3d* box_diagonal = lor_box_offset(lor, universe->cutoff); // the diagonal
	// of the box we need to operate in to catch the entire lor within the cutoff
	// The diagonal will be parralell to the LOR
	vec3d* corner1 = vec_add(lor->center, box_diagonal);
	vec3d* corner2 = vec_sub(lor->center, box_diagonal);
	vec3d* low_corner = three_vec(0,0,0);
	vec3d* high_corner = three_vec(0,0,0);
	low_high_corner(corner1, corner2, low_corner, high_corner);
	free(box_diagonal);
	free(corner1);
	free(corner2);
	vec3d* low_dbl_int = space_to_vox(low_corner, universe);
	vec3d* high_dbl_int = space_to_vox(high_corner, universe);
	free(low_corner);
	free(high_corner);
	int low_x = floor(low_dbl_int->x);
	int low_y = floor(low_dbl_int->y);
	int low_z = floor(low_dbl_int->z);
	int high_x = ceil(high_dbl_int->x);
	int high_y = ceil(high_dbl_int->y);
	int high_z = ceil(high_dbl_int->z);
	free(low_dbl_int);
	free(high_dbl_int);
	if (low_x < 0) {
		low_x = 0;
	}
	if (low_y < 0) {
		low_y = 0;
	}
	if (low_z < 0) {
		low_z = 0;
	}
	if (high_x >= universe->dimensions[0]) {
		high_x = universe->dimensions[0] - 1;
	}
	if (high_y >= universe->dimensions[1]) {
		high_y = universe->dimensions[1] - 1;
	}
	if (high_z >= universe->dimensions[2]) {
		high_z = universe->dimensions[2] - 1;
	}
	for (int i = low_x; i <= high_x; i++) {
		for (int j = low_y; j <= high_y; j++) {
			for (int k = low_z; k <= high_z; k++) {
				// iteration over the entire space in which the lor exists
				vec3d* cur_vox = three_vec(i,j,k);
				vec3d* cur_space = vox_to_space(cur_vox, universe);
				double longitudinal = pt_to_lor_long(lor, cur_space);
				double transverse = pt_to_lor_trans(lor, cur_space);
				free(cur_vox);
				free(cur_space);
				if ((longitudinal < (universe->cutoff * lor->long_uncert))
					&& (transverse < (universe->cutoff * lor->transverse_uncert))) {
					double lon_deviation = longitudinal / lor->long_uncert;
					double trans_deviation = transverse / lor->transverse_uncert;
					if ((universe->cutoff * universe->cutoff) > 
						(lon_deviation * lon_deviation) + (trans_deviation * trans_deviation)) {
						// we are within the processing column (area of useful adding values)
						double lon_normal = long_func(lon_deviation);
						double trans_normal = centered_normal(trans_deviation);
						double total_value = lon_normal * trans_normal;

						int index[3] = {i,j,k};

						pthread_mutex_lock(&volume_lock);
						universe->combiner(universe, total_value, index, attenuation);
						pthread_mutex_unlock(&volume_lock);
					}
				}
			}
		}
	}
}

/*
 * add_lor_plane
 * takes a universe and line of response. Adds the line of response to the
 * universe. This version makes a new iteration plane for
 * each iteration of the third axis
 */
void add_lor_plane(render* universe, lor* lor) {

	int did_we_add_this = 0;
	int did_we_iterate = 0;

	if (universe == NULL || lor == NULL) {
		return;
	}
	// print_lor(stderr, lor);
	// printf("\n");
	// first define the volume in which we will be operating
	vec3d* box_diagonal = lor_box_offset(lor, universe->cutoff); // the diagonal
	// of the box we need to operate in to catch the entire lor within the cutoff
	// The diagonal will be parralell to the LOR
	vec3d* corner1 = vec_add(lor->center, box_diagonal);
	vec3d* corner2 = vec_sub(lor->center, box_diagonal);
	vec3d* low_corner = three_vec(0,0,0);
	vec3d* high_corner = three_vec(0,0,0);
	low_high_corner(corner1, corner2, low_corner, high_corner);
	// creates a corner of low values and a corner of high values
	free(box_diagonal);
	free(corner1);
	free(corner2);
	vec3d* low_dbl_int = space_to_vox(low_corner, universe);
	// the low corner in voxels
	vec3d* high_dbl_int = space_to_vox(high_corner, universe);
	// the high corner in voxels
	free(low_corner);
	free(high_corner);
	int low_x = floor(low_dbl_int->x);
	int low_y = floor(low_dbl_int->y);
	int low_z = floor(low_dbl_int->z);
	int high_x = ceil(high_dbl_int->x);
	int high_y = ceil(high_dbl_int->y);
	int high_z = ceil(high_dbl_int->z);
	free(low_dbl_int);
	free(high_dbl_int);
	if (low_x < 0) {
		low_x = 0;
	}
	if (low_y < 0) {
		low_y = 0;
	}
	if (low_z < 0) {
		low_z = 0;
	}
	if (high_x >= universe->dimensions[0]) {
		high_x = universe->dimensions[0] - 1;
	}
	if (high_y >= universe->dimensions[1]) {
		high_y = universe->dimensions[1] - 1;
	}
	if (high_z >= universe->dimensions[2]) {
		high_z = universe->dimensions[2] - 1;
	}

	// low_x = 0;
	// high_x = universe->dimensions[0] - 1;

	// get the attenuation correction value
	double attenuation = atten_correction(lor);

	// // cos correction for angle to the plane
	// THIS DOES NOT WORK IN A GOOD WAY
	// vec3d z_hat;
	// z_hat.x = 0;
	// z_hat.y = 0;
	// z_hat.z = 1;
	// double correction_angle = vec_angle(&z_hat, lor->dir);
	// attenuation  = attenuation * sin(correction_angle);

	// high and low of each var now defines the maximum extent box, now we
	// define the size of the elipse made by a cylindrical section of the LOR.
	// the section is an elipse with semi-minor axis equal to the LOR radius and
	// a semi-major axis of radius/sin(angle relative to the x_perpedicular plane)

	vec3d x_hat;
	x_hat.x = 1;
	x_hat.y = 0;
	x_hat.z = 0;
	double lor_angle = vec_angle(&x_hat, lor->dir);
	// angle between the xhat vector and the lor
	
	double semi_minor_len = lor->transverse_uncert * universe->cutoff;
	// length of the semi minor axis of the cylinder section (the radius of the
	// LOR cylinder)
	double semi_major_len = (semi_minor_len) / fabs(cos(lor_angle));
	// the semi major axis of the cylinder section of the LOR cylinder

	vec3d* semi_major_dir = three_vec(0.0, lor->dir->y, lor->dir->z);
	// the direction of the semi major axis of the elipse (in the direction of
	// the lor flattened to the yz plane)
	vec3d* semi_major_dir_norm = vec_norm(semi_major_dir);
	// normalized vector pointing in the semi major direction
	free(semi_major_dir);
	vec3d* semi_major = vec_scaler(semi_major_dir_norm, semi_major_len);
	free(semi_major_dir_norm);
	// the semi major axis with direction and amplitude

	vec3d* semi_minor_dir = vec_cross(semi_major, &x_hat); // cross product of
	// the semi major axis and x hat to give a vector perpendicular to the
	// semi major axis
	vec3d* semi_minor_dir_norm = vec_norm(semi_minor_dir);
	free(semi_minor_dir);
	// normalize the semi minor direction
	vec3d* semi_minor = vec_scaler(semi_minor_dir_norm, semi_minor_len);
	free(semi_minor_dir_norm);
	// semi minor axis with direction and magnitude

	vec3d *elipse_box[4];
	elipse_box[0] = vec_add(semi_major, semi_minor);
	elipse_box[1] = vec_sub(semi_major, semi_minor);
	vec3d* negative_semi_major = vec_sub(NULL, semi_major);
	elipse_box[2] = vec_add(negative_semi_major, semi_minor);
	elipse_box[3] = vec_sub(negative_semi_major, semi_minor);

	free(semi_major);
	free(semi_minor);
	free(negative_semi_major);

	// make two corners, they should enclose the entire working area
	vec3d upper_corner;
	upper_corner.y = 0.0;
	upper_corner.z = 0.0;

	for (int i = 0; i < 4; i++) {
		if (elipse_box[i]->y > upper_corner.y) {
			upper_corner.y = elipse_box[i]->y;
		}
		if (elipse_box[i]->z > upper_corner.z) {
			upper_corner.z = elipse_box[i]->z;
		}
		free(elipse_box[i]);
	}

	double upper_corner_y_vox = upper_corner.y * universe->conversion->y;
	double upper_corner_z_vox = upper_corner.z * universe->conversion->z;

	vec3d* high_corner_vox = three_vec(0, upper_corner_y_vox, upper_corner_z_vox);
	vec3d* low_corner_vox = three_vec(0, -upper_corner_y_vox, -upper_corner_z_vox);

	vec_ceil(high_corner_vox);
	vec_floor(low_corner_vox);

	for (int i = low_x; i <= high_x; i++) {
		// find the area that we will be working in. This is the box defined by
		// high and low corner vox centered on the lor at the current x. This
		// volume is then further constrained by the bounding box.
		double x_dist = (i / universe->conversion->x) + universe->least_corner->x - lor->center->x;
		double displacement = x_dist / lor->dir->x; // need to make sure it is not a divide by 0
		vec3d* travel = vec_scaler(lor->dir, displacement);
		vec3d* center = vec_add(travel, lor->center);
		vec3d* center_vox = space_to_vox(center, universe);
		vec_ceil(center_vox);
		int low_elipse_y = (int)low_corner_vox->y + (int)center_vox->y;
		int high_elipse_y = (int)high_corner_vox->y + (int)center_vox->y;
		int low_elipse_z = (int)low_corner_vox->z + (int)center_vox->z;
		int high_elipse_z = (int)high_corner_vox->z + (int)center_vox->z;

		free(travel);
		free(center);

		int y_start;
		int y_end;
		int z_start;
		int z_end;

		if (low_elipse_y < low_y) {
			y_start = low_y;
		} else {
			y_start = low_elipse_y;
		}
		if (low_elipse_z < low_z) {
			z_start = low_z;
		} else {
			z_start = low_elipse_z;
		}
		if (high_elipse_y > high_y) {
			y_end = high_y;
		} else {
			y_end = high_elipse_y;
		}
		if (high_elipse_z > high_z) {
			z_end = high_z;
		} else {
			z_end = high_elipse_z;
		}

		free(center_vox);
		
		for (int j = y_start; j <= y_end; j++) {
			for (int k = z_start; k <= z_end; k++) {
				// iteration over the entire space in which the lor exists
				vec3d* cur_vox = three_vec(i,j,k);
				vec3d* cur_space = vox_to_space(cur_vox, universe);
				double transverse = pt_to_lor_trans(lor, cur_space);
				// transverse distance from the lor to the point
				if (transverse < (universe->cutoff * lor->transverse_uncert)) {
					double longitudinal = pt_to_lor_long(lor, cur_space); 
					// longitudinal distance from the lor center to the point
					if (longitudinal < (universe->cutoff * lor->long_uncert)) {
						// we are within the processing column (area of useful adding values)
						double lon_normal = 1.0;
						if (!long_binary) {
							double long_deviation = longitudinal / lor->long_uncert;
							lon_normal = centered_normal(long_deviation);
						}
						double trans_normal = 1.0;
						if (!trans_binary) {
							double trans_deviation = transverse / lor->transverse_uncert;
							trans_normal = centered_normal(trans_deviation);
						}
						double total_value = lon_normal * trans_normal;

						int index[3] = {i,j,k};
						did_we_add_this = 1;

						pthread_mutex_lock(&volume_lock);
						universe->combiner(universe, total_value, index, attenuation);
						pthread_mutex_unlock(&volume_lock);
					}
				}
				free(cur_vox);
				free(cur_space);
				did_we_iterate = 1;
			}
		}
	}
	free(high_corner_vox);
	free(low_corner_vox);
	if (did_we_iterate) {
		if (!did_we_add_this) {
			pthread_mutex_lock(&volume_lock);
			not_added++;
			iterated++;
			pthread_mutex_unlock(&volume_lock);
		} else {
			pthread_mutex_lock(&volume_lock);
			iterated++;
			pthread_mutex_unlock(&volume_lock);
		}
	}
}

void print_definition(render* rend) {
	printf("Render:\n\tDimensions:\n\t\tx: %i, y: %i, z: %i\n", rend->dimensions[0], rend->dimensions[1], rend->dimensions[2]);
	printf("\tCorners:\n\t\t");
	vec_print(rend->least_corner, stdout);
	printf("\n\t\t");
	vec_print(rend->max_corner, stdout);
	printf("\n");
	printf("\tCombination method: ");
	if (rend->combiner == add_double) {
		printf("addition\n");
	} else if (rend->combiner == add_mult_log) {
		printf("addition with log(multiplication)\n");
	} else if (rend->combiner == add_mult) {
		printf("addition with multiplication\n");
	} else if (rend->combiner == add_binary) {
		printf("addition with binary lor\n");
	} else if (rend->combiner == add_log) {
		printf("addition of log values\n");
	} else {
		printf("unknown pointer\n");
	}
}


void print_volume(FILE* output, render* universe) {
	for (int i = 0; i < universe->dimensions[0]; i++) {
		for (int j = 0; j < universe->dimensions[1]; j++) {
			for (int k = 0; k < universe->dimensions[2]; k++) {
				fprintf(output, "%lf ", universe->volume[i * universe->dimensions[1] * universe->dimensions[2]
															+ j * universe->dimensions[2]
															+ k]);
				
			}
			fprintf(output, "\n");
		}
		// fprintf(output, "\n");
	}
}


void* wrapper_add_lor(void *a) {
	struct _add_lor_union *open = (struct _add_lor_union *)a;
	add_lor(open->universe, open->lor);
	free_lor(open->lor);
	free(open);
	return 0;
}

void* wrapper_add_lor_plane(void *a) {
	struct _add_lor_union *open = (struct _add_lor_union *)a;
	add_lor_plane(open->universe, open->lor);
	free_lor(open->lor);
	free(open);
	return 0;
}


// takes in a universe, location of a lor, location of the update flag, and
// the location of the thread's lock.
// runs as a long term loop getting new lors put into the location it knows to
// look. If the die flag is set to true it stops trying to run and exits after
// running the current lor.
void* long_life_thread_add_lor_plane(void *a) {
	struct _add_lor_union *home = (struct _add_lor_union *)a;
	int die = 0;
	lor* tick_lor = NULL;

	while (!die) {
		// actions to take in the main tick loop
		// check if there is a new LOR to act on
		pthread_mutex_lock(&(home->lock));
		die = home->die;
		if (!home->flag) {
			// nothing this check
			// time to wait for an update
			pthread_cond_wait(&(home->cond), &(home->lock));
			pthread_mutex_unlock(&(home->lock));
			pthread_cond_signal(&main_hold);
		} else {
			// we have a new lor to check out
			tick_lor = home->lor;
			home->flag = 0;
			pthread_mutex_unlock(&(home->lock));
			pthread_cond_signal(&main_hold);
			// now lets run using this new lor
			if (tick_lor != NULL) {
				add_lor_plane(home->universe, tick_lor);
				free_lor(tick_lor);
			}
		}
	}
	return 0;
}





int main(int argc, char const *argv[])
{
	// expected inputs: lor file name, output file name, rendering def file name
	if (argc < 4) {
		if (!strncmp(argv[1], "-h",2) || !strncmp(argv[1], "-H",2)) {
			printf("Kepler's lor renderer: \nThis is designed to be used with");
			printf(" the reverse kinematics code running on a TOPAS simulation.");
			printf("\nThe code expects three files: a lor file, a definitions file");
			printf(" name, and an output file name.\n");
			printf("The lor file should have one line of response for each history.");
			printf(" Each lor should have the struture:\n\thistory_number,");
			printf(" center_x, center_y, center_z, direction_x, dir_y, dir_z,");
			printf(" sigma_longitudianal, sigma_transverse\n\n");
			printf("The output file is just the name you wish for it to have.\n");
			printf("The defintion file is a set of values:\n\tcorner_1_x, ");
			printf("corner_1_y, corner_1_z, corner_2_x, corner_2_y, corner_2_z\n\t");
			printf("voxels_in_x, voxels_in_y, voxels_in_z\n");
			printf("\tLOR_sigma_cutoff\n");
			// add information on selecting combination method and output type.
		}
		printf("Expected 3 inputs, a lor file name, a rendering defintion file name and");
		printf(" an output file name. Use argument -h for more help\n");
		return 1;
	}
	FILE* definition = fopen(argv[2], "r");
	if (definition == NULL) {
		fprintf(stderr, "Unable to open definition file.\n");
		return 1;
	}
	master_copy = read_render_def(definition);
	print_definition(master_copy);
	FILE* input_lor = fopen(argv[1], "r");
	if (input_lor == NULL) {
		fprintf(stderr, "Unable to open LOR file.\n");
		return 1;
	}

	if (argc > 4) {
		for (int i = 1; i < argc; i++) {
			if (!strncasecmp(argv[i], "-g", 2)) {
				if (argc > (i + 1)) {
					// load the geometry
					FILE* in_geo = fopen(argv[i+1], "r");
					if (in_geo != NULL) {
						objects = read_geometry(in_geo);
						fclose(in_geo);
					}
					if (objects == NULL) {
						fprintf(stderr, "WARN: unable to load geometry\n");
					} else {
						geometry_full_tests();
						print_geometry(stdout, objects);
					}
				} else {
					printf("put -g flag before path to geometry file\n");
				}
			}
			// if (!strncasecmp(argv[i], "-lb", 3)) {
			// 	long_func = centered_binary;
			// 	printf("using centered binary for longitudinal\n");
			// }
			// if (!strncasecmp(argv[i], "-ll", 3)) {
			// 	long_func = centered_laplace;
			// 	printf("using centered laplace for longitudinal\n");
			// }
			if (!strncasecmp(argv[i], "-la", 3)) {
				if (argc > (i + 1)) {
					// now set the longitudinal adjustment value
					long_adjust = strtod(argv[i+1], NULL);
					printf("longitudinal adjustment set to %lf\n", long_adjust);
				} else {
					fprintf(stderr, "WARN: -la flag not followed by float\n");
				}
			}
			if (!strncasecmp(argv[i], "-ta", 3)) {
				if (argc > (i + 1)) {
					// now set the transverse adjustment value
					trans_adjust = strtod(argv[i+1], NULL);
					printf("transverse adjustment set to %lf\n", trans_adjust);
				} else {
					fprintf(stderr, "WARN: -ta flag not followed by float\n");
				}
			}
			if (!strncasecmp(argv[i], "-lb", 3)) {
				// now set the longitudinal to binary
				long_binary = 1;
				printf("longitudinal set to binary\n");
			}
			if (!strncasecmp(argv[i], "-tb", 3)) {
				// now set the transverse to binary
				trans_binary = 1;
				printf("transverse set to binary\n");			}
		}
	}
	if (long_func == NULL) {
		long_func = centered_normal;
	}

	if (pthread_mutex_init(&volume_lock, NULL)) {
		fprintf(stderr, "pthread: unable to make volume lock, exiting.\n");
		return(1);
	}

	for (int i = 0; i < MAX_THREADS; i++) {
		// set the array of threads to be empty
		tid[i] = -1;
		working_lors[i].lor = NULL;
		working_lors[i].die = 0;
		working_lors[i].flag = 0;
		working_lors[i].universe = master_copy;
		pthread_cond_init(&(working_lors[i].cond), NULL);
		if (pthread_mutex_init(&(working_lors[i].lock), NULL)) {
			fprintf(stderr, "pthread: unable to make update lock, exiting.\n");
			return(1);
		}
		pthread_create(&tid[i], NULL, long_life_thread_add_lor_plane, &(working_lors[i]));
	}
	pthread_cond_init(&main_hold, NULL);
	pthread_mutex_init(&main_lock, NULL);

	uint iteration = 0;
	lor* operative_lor = read_lor(input_lor);

	int cur_thread = 0;
	while (operative_lor != NULL) {
		int cont = 0;
		while (!cont) {
			for (int i = 0; (i < MAX_THREADS) && (cont == 0); i++) {
				// work our way through the threads, updating the next avaliable one
				// with our current LOR.
				pthread_mutex_lock(&(working_lors[i].lock));
				// lock it so we can have a look
				if (!working_lors[i].flag) {
					// the thread already grabbed a lor to work on, time to point to
					// the new one
					working_lors[i].lor = operative_lor;
					working_lors[i].flag = 1;
					// wake the thread back up if it had been suspended
					pthread_cond_signal(&(working_lors[i].cond));
					cont = 1;
				}
				pthread_mutex_unlock(&(working_lors[i].lock));
			}
			if (!cont) {
				pthread_mutex_lock(&main_lock);
				pthread_cond_wait(&main_hold, &main_lock);
				pthread_mutex_unlock(&main_lock);
			}
		}
		
		// if (tid[cur_thread] != -1) {
		// 	pthread_join(tid[cur_thread], NULL);
		// 	tid[cur_thread] = -1;
		// }
		// // add_lor(master_copy, operative_lor);
		// add_lor_plane(master_copy, operative_lor);

		// lets make this multithreaded
		// struct _add_lor_union *arguments = (struct _add_lor_union *)malloc(sizeof(struct _add_lor_union));
		// arguments->universe = master_copy;
		// arguments->lor = operative_lor;
		// pthread_create(&tid[cur_thread], NULL, wrapper_add_lor_plane, arguments);
		// pthread_create(&tid[cur_thread], NULL, wrapper_add_lor_plane, arguments);

		operative_lor = read_lor(input_lor);
		iteration++;
		if (cur_thread + 1 < MAX_THREADS) {
			cur_thread++;
		} else {
			cur_thread = 0;
		}
		if (100000 * (iteration / 100000) == iteration) {
			printf("iteration %i\n", iteration);
		}
	}

	FILE* output = fopen(argv[3], "w");

	print_volume(output, master_copy);

	printf("Gaussian lookup table:\n\tFails: %llu\n\tSuccesses: %llu\n", gauss_lookup_fail, gauss_lookup_success);
	printf("LORs not added:\n\t%li, fractionally %lf\n", not_added, (double)not_added / (double)iterated);

	pthread_mutex_destroy(&volume_lock);

	return 0;
}
