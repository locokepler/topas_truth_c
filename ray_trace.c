#include "vector_ops.h"
#include "llist.h"
#include "ray_trace.h"
#include <stdlib.h>


/* 
 * Simple ray tracing code to allow for attenuation correction in the LOR
 * rendering operations.
 * Allows for the building of a combined shape using subshapes. The current
 * allowed shapes are the sphere, cylinder, and rectangular prism. The shapes
 * are defined by their centerpoint, their dimentions, and in the case of the
 * cylinder: what axis the circle cross-section is normal to. Also includes
 * information about the attenuation coefficent, however it will produce
 * undefined behavior if the shape has different attenuation coefficents at
 * different locations.
 * 
 * Allows for overlapping geometry, as the method of ray tracing is taking a
 * ray, checking if it is inside the geometry, then calculating the next
 * intersection point in the direction of travel. When it exits the geometry
 * object it is currently in, it checks if it is inside of any other objects.
 * If it is, then the distance to the next exit is found and added to the
 * running total.
 * 
 * If a ray does not intersect a shape in the given direction then a NULL
 * pointer is returned by the intersection function of the shape.
 * 
 * 
 * 
 */

// builds a ray structure
ray* ray_build(vec3d* pos, vec3d* dir) {
    if ((pos == NULL) || (dir = NULL)) {
        return NULL;
    }
    ray* new = (ray*)malloc(sizeof(ray));
    if (new == NULL) {
        fprintf(stderr, "make_ray: malloc failed, likely fatal\n");
        return NULL;
    }
    new->pos = pos;
    new->dir = dir;
    return new;
}

// copy a ray structure
ray* ray_copy(ray* src) {
    if (src == NULL) {
        return NULL;
    }
    vec3d* pos_cpy = vec_copy(src->pos);
    vec3d* dir_cpy = vec_copy(src->dir);
    return ray_build(pos_cpy, dir_cpy);
}

// frees a ray structure
void ray_free(ray* src) {
    if (src == NULL) {
        return;
    }
    if (src->dir != NULL) {
        free(src->dir);
    }
    if (src->pos != NULL) {
        free(src->pos);
    }
}

traversal* traversal_build(vec3d* intersect, double t) {
    if (intersect == NULL) {
        return NULL;
    }
    traversal* new = (traversal*)malloc(sizeof(traversal));
    if (new == NULL) {
        fprintf(stderr, "traversal_build: malloc failed, likely fatal\n");
        return NULL;
    }
    new->intersection = intersect;
    new->t = t;
    return new;
}

traversal* traversal_copy(traversal* src) {
    return traversal_build(src->intersection, src->t);
}

void traversal_free(traversal* src) {
    free(src->intersection);
    free(src);
}


/*
 * Takes the vector and adjusts the values to transfer them to the given axis.
 * The change is moving the z axis to the given location. 
 * (1 = x, 2 = y, 3=z(do nothing)). Negative values are also accepted. A
 * negative value is interpreted as a rotation in the opposite direction. If it 
 * was able to do the transformation a 1 is returned. If it failed in any way it
 * returns a 0.
 */
int coord_transfer(vec3d* coord, int axis) {
    if ((coord == NULL) || (axis > 3) || (axis < -3) || (axis == 0)) {
        return 0;
    }
    if ((axis == 1) || (axis == -1)) {
        // move the z axis onto the x axis 
        double x_old = coord->x;
        double z_old = coord->z;
        coord->x = ((double)axis) * z_old;
        coord->z = -((double)axis) * x_old;
        return 1;
    } else if ((axis == 2) || (axis == -2)) {
        // move the z axis onto the y axis
        double y_old = coord->y;
        double z_old = coord->z;
        coord->y = ((double)axis) * z_old;
        coord->z = -((double)axis) * y_old;
        return 1;
    } else {
        // axis is 3, for 3 we keep the z where it is, same for -3
    }
    return 1;
}

/*
 * finds if the given ray (propagate) intersects a given plane. The plane is
 * defined by the normal axis (axis = 0,1,2 = x,y,z normal) and the size in
 * floats. The normal axis is to be treated as the z axis, with the size 
 * defining the size along the x (size[0]) and y (size[1]) directions. This plane
 * is centered at the location given by center.
 */
traversal* plane_intersect_rec(float size[2], int axis, vec3d* center_src, ray* propagate_src) {
    if ((size == NULL) || (center_src == NULL) || (propagate_src == NULL)) {
        return NULL;
    }
    // rotate everything so that it is all lined up
    vec3d* center = vec_copy(center);
    ray* propagate = ray_copy(propagate_src);
    coord_transfer(center, axis);
    coord_transfer(propagate->pos, axis);
    coord_transfer(propagate->dir, axis);
    // find when the intersection should happen (when ray crosses z of plane)
    double t = (center->z - propagate->pos->z) / propagate->dir->z;
    double x_cross = propagate->pos->x + (propagate->dir->x * t);
    double x_hl = 0.5 * size[0];
    if ((x_cross > (x_hl + center->x)) || (x_cross < (-x_hl + center->x))) {
        // we missed the x edge of the plane, no intersection
        return NULL;
    }
    double y_cross = propagate->pos->y + (propagate->dir->y * t);
    double y_hl = 0.5 * size[1];
    if ((y_cross > (y_hl + center->y)) || (y_cross < (-y_hl + center->y))) {
        return NULL;
    }
    vec3d* intersection = three_vec(x_cross, y_cross, center->z);
    // now that we have the intersection we need to rotate everything back to
    // its original position (negative axis value)
    coord_transfer(intersection, -axis);
    free(center);
    ray_free(propagate);
    traversal* output = traversal_build(intersection, t);
    return output;
}

/*
 * finds if the given ray (propagate) intersects a given plane. The plane is
 * defined by the normal axis (axis = 0,1,2 = x,y,z normal) and the size in
 * floats. The normal axis is to be treated as the z axis, with the size 
 * defining the radius as size. This plane
 * is centered at the location given by center.
 */

traversal* plane_intersect_circle(float size, int axis, vec3d* center_src, ray* propagate_src) {
    if ((center_src == NULL) || (propagate_src == NULL)) {
        return NULL;
    }
    // rotate everything so that it is all lined up
    vec3d* center = vec_copy(center);
    ray* propagate = ray_copy(propagate_src);
    coord_transfer(center, axis);
    coord_transfer(propagate->pos, axis);
    coord_transfer(propagate->dir, axis);
    // find when the intersection should happen (when ray crosses z of plane)
    double t = (center->z - propagate->pos->z) / propagate->dir->z;
    double x_cross = propagate->pos->x + (propagate->dir->x * t);
    double y_cross = propagate->pos->y + (propagate->dir->y * t);
    double rad2 = size * size;
    double x_disp = x_cross - center->x;
    double y_disp = y_cross - center->y;
    double disp2 = (x_disp * x_disp) + (y_disp * y_disp);
    if (rad2 < disp2) {
        // intersection not within the given radius
        return NULL;
    }
    vec3d* intersection = three_vec(x_cross, y_cross, center->z);
    // now that we have the intersection we need to rotate everything back to
    // its original position (negative axis value)
    coord_transfer(intersection, -axis);
    free(center);
    ray_free(propagate);
    traversal* output = traversal_build(intersection, t);
    return output;
}

traversal* exit_rectangular_prism(ray* path, shape* prism) {
    if ((path == NULL) || (prism == NULL) || (prism->type != REC_PRISM)) {
        return NULL;
    }
    double long_path = -1.0;
    double short_path = -1.0;
    // check the distance for all six faces.
    for (int i = 0; i < 3; i++) {
        // check the face in the + axis direction and the -axis direction
        // each face has the dimensions of the two axis not using i
        float size[2];
        if (i < 2) {
            size[0] = prism->dimentions[i + 1];
        } else {
            size[0] = prism->dimentions[0];
        }
        if (i < 1) {
            size[1] = prism->dimentions[i + 2];
        } else {
            size[1] = prism->dimentions[i - 1];
        }
        // get the offset of the plane we will check from the center of the shape
        double offset = size[i] * 0.5;
        // centerpoints of the planes
        vec3d* direction_p = three_vec(prism->position[0] + (offset * (i == 0)),
                                        prism->position[1] + (offset * (i == 1)),
                                        prism->position[2] + (offset * (i == 2)));
        vec3d* direction_m = three_vec(prism->position[0] - (offset * (i == 0)),
                                        prism->position[1] - (offset * (i == 1)),
                                        prism->position[2] - (offset * (i == 2)));
        traversal* exit1 = plane_intersect_rec(size, i + 1, direction_p, path);
        traversal* exit2 = plane_intersect_rec(size, i + 1, direction_m, path);
        // now check if either has an intersection at all,
        if ((exit1 != NULL) && (exit1->t > 0)) {
            if (long_path < 0) {
                long_path = exit1->t;
            } else {
                if (exit1->t > long_path) {
                    short_path = long_path;
                    long_path = exit1->t;
                } else {
                    short_path = exit1->t;
                }
            }
        }
        if ((exit2 != NULL) && (exit2->t > 0)) {
            if (long_path < 0) {
                long_path = exit2->t;
            } else {
                if (exit2->t > long_path) {
                    short_path = long_path;
                    long_path = exit2->t;
                } else {
                    short_path = exit2->t;
                }
            }
        }
    }
}