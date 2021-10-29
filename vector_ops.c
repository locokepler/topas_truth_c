#include <math.h>
#include <stdlib.h>
#include "vector_ops.h"

// defines a new 3-vector
vec3* three_vec(double x, double y, double z) {
	vec3* new = (vec3*)malloc(sizeof(vec3));
	new->x = x;
	new->y = y;
	new->z = z;
	return new;
}

// determines the distance between two vectors
double vec_dist(vec3* a, vec3* b) {
	vec3* new = vec_sub(a, b);
	double distance = vec_mag(new);
	// free the new vector created to measure the magnitude
	free(new);
	return distance;
}

// returns the magnitude of a 3-vector
double vec_mag(vec3* vector) {
	if (vector == NULL) {
		return NAN;
	}
	double x2 = vector->x * vector->x;
	double y2 = vector->y * vector->y;
	double z2 = vector->z * vector->z;
	return sqrt(x2 + y2 + z2);
}

// returns the dot product of two 3-vectors
double vec_dot(vec3* a, vec3* b) {
	if (a == NULL || b == NULL) {
		return NAN;
	}
	double first = a->x * b->x;
	double second = a->y * b->y;
	double third = a->z * b->z;
	return first + second + third;
}

// returns the addition of two 3-vectors as a new 3-vector
vec3 *vec_add(vec3* a, vec3* b) {
	if (a == NULL || b == NULL) {
		return NULL;
	}
	double x = a->x + b->x;
	double y = a->y + b->y;
	double z = a->z + b->z;
	return three_vec(x, y, z);
}

// returns the subtraction of the second vector from the first
vec3 *vec_sub(vec3* a, vec3* b) {
	if (a == NULL || b == NULL) {
		return NULL;
	}
	double x = a->x - b->x;
	double y = a->y - b->y;
	double z = a->z - b->z;
	return three_vec(x, y, z);
}

// returns the cross product of the two vectors a X b
vec3 *vec_cross(vec3* a, vec3* b) {
	if (a == NULL || b == NULL) {
		return NULL;
	}
	double x = (a->y * b->z) - (a->z * b->y);
	double y = (a->x * b->z) - (a->z * b->x);
	double z = (a->x * b->y) - (a->y * b->x);
	return three_vec(x,y,z);
}

// calculates the angle between two 3 vectors
double vec_angle(vec3* a, vec3* b) {
	if (a == NULL || b == NULL) {
		return -1.;
	}
	// arccos((a * b) / (||a|| ||b||)) = theta
	double a_dot_b = vec_dot(a, b);
	double norms = vec_mag(a) * (vec_mag(b));
	return acos(a_dot_b / norms);
}

// makes a new copy of a vector
vec3* vec_copy(vec3* a) {
	if (a == NULL) {
		return NULL;
	}
	return three_vec(a->x, a->y, a->z);
}