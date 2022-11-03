#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"
#include <pthread.h>
#include "ray_trace.h"
#include "llist.h"

#define MAX_THREAD_CALLS 4

#define PHI_LORS

pthread_t tid[MAX_THREAD_CALLS];
pthread_mutex_t volume_lock;

geometry* objects = NULL;

double (*long_func)(double, double) = NULL;

double long_adjust = 1.0;
double trans_adjust = 1.0;

typedef struct _intvec {
	int a;
	int b;
	int c;
} int_vec;

/* Takes a file of lines of response and builds them into a voxelized image
 * by taking the LOR and stepping through each voxel it crosses. The LOR weight
 * at the midpoint of the LOR's journey through the voxel is added (times an
 * attenuation weight). This lacks a transverse resolution behavior.
 * 
 * Takes the same definition and .lor files as lor_render1
 */


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
void add_double(render* u, double x, int* indexes, double arg) {
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += x;
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
double centered_normal(double sigma, double x) {
	double fraction = x / sigma;
	double exponent = -0.5 * (fraction * fraction);
	return exp(exponent);
}
double centered_binary(double sigma, double x) {
	return 1.0;
}

double centered_laplace(double sigma, double x) {
	double fraction = x / sigma;
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

int add_along_ray(render* universe, ray* travel, double atten, double max_dist, double long_uncert) {
	if ((universe == NULL) || (travel == NULL)) {
		return 0;
	}

    // first we need to know what voxel we are in
    vec3d* start_vox = space_to_vox(travel->pos, universe); // our location in voxels
    vec_floor(start_vox); // move to along index values
    vec3d* box_start_space = vox_to_space(start_vox, universe); // location of our
	// voxel location in space
    vec3d* voxel_size_cm = three_vec(1/universe->conversion->x,
									1/universe->conversion->y,
									1/universe->conversion->z); // dimentions of a voxel in cm

    float voxel_dims[3] = {voxel_size_cm->x, voxel_size_cm->y, voxel_size_cm->z}; // the size of a voxel (in cm)
    float voxel_loc[3] = {box_start_space->x + 0.5 * voxel_size_cm->x,
                            box_start_space->y + 0.5 * voxel_size_cm->y,
                            box_start_space->z + 0.5 * voxel_size_cm->z};
							// voxel location with voxel centered between the indices

    shape* voxel = shape_build(REC_PRISM, voxel_loc, voxel_dims, 0, 1.0);
	// build the shape of the voxel

	int index[3] = {floor(start_vox->x), floor(start_vox->y), floor(start_vox->z)};
	double total_dist = 0;

	while ((total_dist < max_dist)
			 && ((index[0] >= 0) && (index[0] < universe->dimensions[0]))
			 && ((index[1] >= 0) && (index[1] < universe->dimensions[1]))
			 && ((index[2] >= 0) && (index[2] < universe->dimensions[2]))) {
		// we must be within the max dist, and also not outside of the correct index
		int cross;
		// how far did we go in the voxel, and where did we exit?
		traversal* exit = exit_rectangular_prism(travel, voxel, &cross);
		if (exit == NULL) {
			fprintf(stderr, "add_along_ray: ERROR! failed to intersect voxel!\n");
			return 0;
		}
		vec3d* midpoint_adj = vec_scaler(travel->dir, -0.5 * exit->t);
		vec3d* midpoint = vec_sub(exit->intersection, midpoint_adj);
		total_dist = vec_dist(midpoint, travel->pos); // also by def the longitudinal dist
		double lon_gauss = long_func(long_uncert, total_dist);
		lon_gauss *= atten;
		// add to the given volume
		pthread_mutex_lock(&volume_lock);
		universe->combiner(universe, lon_gauss, index, atten);
		pthread_mutex_unlock(&volume_lock);
		// find where next to look. Need the new index
		// give the next checked location as being the voxelized location of
		// exit + dt
		double dt = 0.0001; // should be small to not skip voxels
		vec3d* dir_dt = vec_scaler(travel->dir, dt);
		vec3d* checkpoint = vec_add(exit->intersection, dir_dt);
		// difference between checkpoint voxel number and current index
		vec3d* vox_checkpoint = space_to_vox(checkpoint, universe);
		vec_floor(vox_checkpoint);
		int dx = vox_checkpoint->x - index[0];
		int dy = vox_checkpoint->y - index[1];
		int dz = vox_checkpoint->z - index[2];
		if ((dx == 0) && (dy == 0) && (dz == 0)) {
			fprintf(stderr, "add_along_ray: ERROR! failed to move!\n");
			return 0;
		}
		// set the index to the new voxel location
		index[0] = vox_checkpoint->x;
		index[1] = vox_checkpoint->y;
		index[2] = vox_checkpoint->z;
		// shift the voxel location by change in index.
		voxel->pos[0] += voxel->dim[0] * dx;
		voxel->pos[1] += voxel->dim[1] * dy;
		voxel->pos[2] += voxel->dim[2] * dz;
		traversal_free(exit);
		free(midpoint_adj);
		free(midpoint);
		free(dir_dt);
		free(checkpoint);
		free(vox_checkpoint);
	}

	free(start_vox);
	free(box_start_space);
	free(voxel_size_cm);
	shape_free(voxel);

	return 1;
}

/*
 * add_lor
 * takes a universe and line of response. Adds the line of response to the
 * universe. It takes a direction and propagates the LOR along the direction
 * until it hits the sigma distance limit. Turns around at the center and
 * repeats in the other direction. For each voxel it adds the weight of the line
 * at the midpoint of its trip through the voxel.
 */
void add_lor(render* universe, lor* lor) {
	if (universe == NULL || lor == NULL) {
		return;
	}
	// get the attenuation correction value
	double attenuation = atten_correction(lor);
	// check if we are in the image volume.
	float dims_of_universe[3] = {(float)(universe->dimensions[0]) / universe->conversion->x,
									(float)(universe->dimensions[1]) / universe->conversion->y,
									(float)(universe->dimensions[2]) / universe->conversion->z};
	vec3d* center_offset_vox = three_vec(((float)(universe->dimensions[0])) * 0.5, 
											((float)(universe->dimensions[1])) * 0.5, 
											((float)(universe->dimensions[2])) * 0.5);
	vec3d* center_offset_space = vox_to_space(center_offset_vox, universe);
	float center_of_universe[3] = {center_offset_space->x,
									center_offset_space->y,
									center_offset_space->z};
	shape* volume = shape_build(REC_PRISM, center_of_universe, dims_of_universe, 0,1);
    vec3d* go = vec_copy(lor->dir); // is automatically normalized
    vec3d* start_space = vec_copy(lor->center);
    ray* travel = ray_build(start_space, go);
	int cross = 0;
	traversal* exit = exit_rectangular_prism(travel, volume, &cross);
	free(center_offset_vox);
	free(center_offset_space);
	shape_free(volume);
	if (exit != NULL) {
		// In the current direction of the LOR we interact with the volume
		// add the lor in this direction. Start at when calculation of this ray's
		// propagation inside the volume started. This is at the center of the
		// lor for a partial crossing, at the entrance to the volume for a full
		// crossing.
		vec3d* displacement = vec_scaler(travel->dir, -exit->t + 0.0001);
		// has a small adjustment (1 um) to keep it inside of the volume
		vec3d* start = vec_add(exit->intersection, displacement);
		ray* path = ray_build(start, vec_copy(travel->dir));
		// add along the current path
		add_along_ray(universe, path, attenuation, lor->long_uncert * universe->cutoff, lor->long_uncert);
		free(displacement);
		ray_free(path);
		if (!cross) {
			// we did not cross the entire volume, so we need to do this again
			// heading the other direction
			vec3d* displacement = vec_scaler(travel->dir, -exit->t - 0.0001);
			// has a small adjustment (1 um) to keep it inside of the volume
			vec3d* start = vec_add(exit->intersection, displacement);
			vec3d* reverse = vec_sub(NULL, travel->dir);
			ray* path = ray_build(start, reverse);
			// add along reverse path
			add_along_ray(universe, path, attenuation, lor->long_uncert * universe->cutoff, lor->long_uncert);
			free(displacement);
			ray_free(path);
			traversal_free(exit);
		}
	} else {
		vec3d* reverse = vec_sub(NULL, travel->dir);
		free(travel->dir);
		travel->dir = reverse;
		exit = exit_rectangular_prism(travel, volume, &cross);
		if (exit != NULL) {
			vec3d* displacement = vec_scaler(travel->dir, -exit->t + 0.0001);
			vec3d* start = vec_add(exit->intersection, displacement);
			ray* path = ray_build(start, vec_copy(travel->dir));
			add_along_ray(universe, path, attenuation, lor->long_uncert * universe->cutoff, lor->long_uncert);
			free(displacement);
			ray_free(path);
			traversal_free(exit);
		} else {
			// this ray apparently does not intersect the volume!
		}

	}
	ray_free(travel);


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

struct _add_lor_union {
	render* universe;
	lor* lor;
};

void* wrapper_add_lor(void *a) {
	struct _add_lor_union *open = (struct _add_lor_union *)a;
	add_lor(open->universe, open->lor);
	free_lor(open->lor);
	free(open);
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
			if (!strncasecmp(argv[i], "-lb", 3)) {
				long_func = centered_binary;
				printf("using centered binary for longitudinal\n");
			}
			if (!strncasecmp(argv[i], "-ll", 3)) {
				long_func = centered_laplace;
				printf("using centered laplace for longitudinal\n");
			}
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
		}
	}
	if (long_func == NULL) {
		long_func = centered_normal;
	}

	if (pthread_mutex_init(&volume_lock, NULL)) {
		fprintf(stderr, "pthread: unable to make volume lock, exiting.\n");
		return(1);
	}

	for (int i = 0; i < MAX_THREAD_CALLS; i++) {
		// set the array of threads to be empty
		tid[i] = -1;
	}

	uint iteration = 0;
	lor* operative_lor = read_lor(input_lor);

	int cur_thread = 0;
	while (operative_lor != NULL) {
		if (tid[cur_thread] != -1) {
			pthread_join(tid[cur_thread], NULL);
			tid[cur_thread] = -1;
		}
		// add_lor(master_copy, operative_lor);

		// // lets make this multithreaded
		// struct _add_lor_union *arguments = (struct _add_lor_union *)malloc(sizeof(struct _add_lor_union));
		// arguments->universe = master_copy;
		// arguments->lor = operative_lor;
		// pthread_create(&tid[cur_thread], NULL, wrapper_add_lor, arguments);
		// pthread_create(&tid[cur_thread], NULL, wrapper_add_lor_plane, arguments);
        add_lor(master_copy, operative_lor);
        free_lor(operative_lor);

		operative_lor = read_lor(input_lor);
		iteration++;
		if (cur_thread + 1 < MAX_THREAD_CALLS) {
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

	pthread_mutex_destroy(&volume_lock);

	return 0;
}
