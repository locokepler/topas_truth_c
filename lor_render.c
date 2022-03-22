#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"
#include "lor.h"

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
	worked = fscanf(input, "%lf", &longitudinal);
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
	vec3* cleanup = three_vec(dir_x, dir_y, dir_z);
	new->dir = vec_norm(cleanup);
	free(cleanup);
	new->long_uncert = longitudinal;
	new->transverse_uncert = transverse;


	return new;
}

// takes a render (for the array), a double to add to it, and the location to
// do so. The indexes array should be of length 3, an x, y, and z coord.
void add_double(render* u, double x, int* indexes) {
	u->volume[indexes[0] * u->dimensions[1] * u->dimensions[2]
				+ indexes[1] * u->dimensions[2]
				+ indexes[2]] += x;
	return;
}


// takes two corner defintions and rebuilds them so that they have one of all
// low values and one of all high values. Makes iteration simpler.
// All operations done in situe with given locations
void low_high_corner(vec3* corner1, vec3* corner2, vec3* low_corner, vec3* high_corner) {
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
	vec3* corner1 = three_vec(x, y, z);
	worked = fscanf(input, "%lf,%lf,%lf", &x, &y, &z); // render cube corner 2
	vec3* corner2 = three_vec(x, y, z);
	vec3* low_corner = three_vec(0,0,0);
	vec3* high_corner = three_vec(0,0,0);
	low_high_corner(corner1, corner2, low_corner, high_corner);
	free(corner1);
	free(corner2);
	int a, b, c;
	worked = fscanf(input, "%i,%i,%i", &a, &b, &c); // gets the number of voxels in each dimension
	
	int cutoff;
	worked = fscanf(input, "%i", &cutoff);

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


	// add stuff to handle various choices of conversion here when needed.

	new->combiner = add_double;
	return new;
}

/* 
 * space_to_vox
 * converts from spacial coordiantes (typically in cm) to the voxel coordinates.
 * This conversion does not do the final change to integeriaztion. The values
 * returns will likely still need to be rounded or similar to define the final
 * values.
 */
vec3* space_to_vox(vec3* a, render* universe) {
	if (a == NULL || universe == NULL) {
		return NULL;
	}
	// first handle the offset
	vec3* offsetless = vec_sub(a, universe->least_corner);
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
vec3* vox_to_space(vec3* a, render* universe) {
	if (a == NULL || universe == NULL) {
		return NULL;
	}
	double x = a->x / universe->conversion->x;
	double y = a->y / universe->conversion->y;
	double z = a->z / universe->conversion->z;

	vec3* multiplied = three_vec(x, y, z);
	vec3* offset = vec_add(multiplied, universe->least_corner);
	free(multiplied);

	return offset;
}

/*
 * pt_to_lor_long
 * Gives the longitudinal distance from center of the LOR to the given point.
 * Done using a simple dot product using the fact the direction vector is
 * a normalized vector
 */
double pt_to_lor_long(lor* line, vec3* pt) {
	if (line == NULL || pt == NULL) {
		return -1;
	}
	vec3* to_point = vec_sub(pt, line->center);
	double dist = vec_dot(to_point, line->dir);
	free(to_point);
	return fabs(dist);
}

/*
 * pt_to_lor_trans
 * Gives the transverse distance from the center of the lor to the given point.
 * Done using a cross product and the fact the direction vector is normalized
 */
double pt_to_lor_trans(lor* line, vec3* pt) {
	if (line == NULL || pt == NULL) {
		return -1;
	}
	vec3* to_point = vec_sub(pt, line->center);
	vec3* cross = vec_cross(to_point, line->dir);
	double dist = vec_mag(cross);
	free(to_point);
	free(cross);
	return fabs(dist);
}

/*
 * lor_box_offset
 * Takes a lor and make the offset for half of the processing box. The point
 * is defined as being at cutoff * sigma_l + (cutoff * sigma_t in the x, y, and
 * z direction that the lor points in). This guarrentees that the entire volume
 * of the lor within the cutoff band is contained in the volume.
 */
vec3* lor_box_offset(lor* lor, double cutoff_sigma) {
	vec3* longitudinal = vec_scaler(lor->dir, lor->long_uncert * cutoff_sigma);
	double x = (lor->dir->x > 0) - (lor->dir->x < 0); // 1 if heading in x dir, -1 if in -x dir 0 if 0
	double y = (lor->dir->y > 0) - (lor->dir->y < 0);
	double z = (lor->dir->z > 0) - (lor->dir->z < 0);
	vec3* trans_dir = three_vec(x, y, z);
	vec3* trans_vec = vec_scaler(trans_dir, cutoff_sigma * lor->transverse_uncert);
	vec3* total = vec_add(longitudinal, trans_vec);
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

void add_lor(render* universe, lor* lor) {
	// first define the volume in which we will be operating
	vec3* box_diagonal = lor_box_offset(lor, universe->cutoff);
	vec3* corner1 = vec_add(lor->center, box_diagonal);
	vec3* corner2 = vec_sub(lor->center, box_diagonal);
	vec3* low_corner = three_vec(0,0,0);
	vec3* high_corner = three_vec(0,0,0);
	low_high_corner(corner1, corner2, low_corner, high_corner);
	free(corner1);
	free(corner2);
	vec3* low_dbl_int = space_to_vox(low_corner, universe);
	vec3* high_dbl_int = space_to_vox(high_corner, universe);
	free(low_corner);
	free(high_corner);
	int low_x = floor(low_dbl_int->x);
	int low_y = floor(low_dbl_int->y);
	int low_z = floor(low_dbl_int->z);
	int high_x = ceil(high_dbl_int->x);
	int high_y = ceil(high_dbl_int->y);
	int high_z = ceil(high_dbl_int->z);
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
				vec3* cur_vox = three_vec(i,j,k);
				vec3* cur_space = vox_to_space(cur_vox, universe);
				double longitudinal = pt_to_lor_long(lor, cur_space);
				double transverse = pt_to_lor_trans(lor, cur_space);
				if ((longitudinal < (universe->cutoff * lor->long_uncert))
					&& (transverse < (universe->cutoff * lor->transverse_uncert))) {
					// we are within the processing column (area of useful adding values)
					double lon_deviation = longitudinal / lor->long_uncert;
					double trans_deviation = transverse / lor->transverse_uncert;
					double lon_normal = centered_normal(lor->long_uncert, lon_deviation);
					double trans_normal = centered_normal(lor->transverse_uncert, trans_deviation);
					double total_value = lon_normal * trans_normal;

					int index[3] = {i,j,k};

					universe->combiner(universe, total_value, index);
					
				}
			}
		}
	}
	
}













int main(int argc, char const *argv[])
{
	// expected inputs: lor file name, output file name, rendering def file name
	if (argc != 4) {
		if (!strncmp(argv[1], "-h",2) || !strncmp(argv[1], "-H",2)) {
			printf("Kepler's lor renderer: \nThis is designed to be used with");
			printf(" the reverse kinematics code running on a TOPAS simulation.");
			printf("\nThe code expects three files: a lor file, an output file ");
			printf(" name, and a definitions file.\n");
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
		printf("Expected 3 inputs, a lor file name, an output file name and");
		printf(" a rendering defintion file name. Use argument -h for more help\n");
		return 1;
	}
	FILE* definition = fopen(argv[2], "r");
	master_copy = read_render_def(definition);
	return 0;
}
