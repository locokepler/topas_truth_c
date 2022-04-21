#ifndef vector_ops_h
#define vector_ops_h

#include <stdio.h>

typedef struct vec3_double {
	double x;
	double y;
	double z;
} vec3d;

// defines a new 3-vector
vec3d* three_vec(double x, double y, double z);

// returns the magnitude of a 3-vector
double vec_mag(vec3d* vector);

// returns the dot product of two 3-vectors
double vec_dot(vec3d* a, vec3d* b);

// returns the addition of two 3-vectors as a new 3-vector
vec3d *vec_add(vec3d *a, vec3d* b);

// returns the subtraction of the second vector from the first
// if the first vector is NULL acts as negating
vec3d *vec_sub(vec3d* a, vec3d* b);

// returns the cross product of the two vectors a X b as a new vector
vec3d *vec_cross(vec3d* a, vec3d* b);

// calculates the angle between two 3 vectors in radians
double vec_angle(vec3d* a, vec3d* b);

// determines the distance between two vectors
double vec_dist(vec3d* a, vec3d* b);

// makes a new copy of a vector
vec3d* vec_copy(vec3d* a);

// prints a vector as the three values
void vec_print(vec3d* a, FILE* output);

// normalizes the given vector, returns as a new vector structure
vec3d* vec_norm(vec3d* a);

// multiplies a vector by a scalar
vec3d* vec_scaler(vec3d* a, double b);

#endif