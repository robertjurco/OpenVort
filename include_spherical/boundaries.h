#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "vec3_math.h"

#define DEG2RAD(X) ((X)*M_PI/180.0)

/*******************************************************************************
 ******************************** DOMAIN BOX ***********************************
 ******************************************************************************/

/*
	Enumberator of boundary surfaces: NOT_A_FACE (-1), INNER, OUTER.
*/
typedef enum _sf {
	NOT_A_SURFACE = -1, // not a surface
	INNER,
	OUTER
} boundary_surfaces;

/*
	Enumberator of node types: EMPTY (-1), FREE, PINNED, PINNED_SLIP.
*/
typedef enum _ns_e {
  EMPTY = -1, //no point here
  FREE, //ordinary vortex point
  PINNED, //pinned on the wall
  PINNED_SLIP //pinned on the wall, but can slip
} node_status_t;

/*
	Status of the node and on which wall it is pinned.
*/
struct node_status {
  node_status_t status; // Node status: EMPTY (-1), FREE, PINNED, PINNED_SLIP.
  boundary_surfaces pin_surface; // On which wall the node is pinned: NOT_A_FACE (-1), INNER, OUTER.
};

/*
	Size of the domain sphere, angular size is in both cases (-angle, +angle).
*/
struct domain {
	double inner_radius;
	double outer_radius;
	double azimut_angle;
	double polar_angle;
};

/*******************************************************************************
 ********************** PERIODIC AND MIRROR GEOMETRIES *************************
 ******************************************************************************/

struct vec3 surface_normal(const struct vec3* where, int surface);

struct segment seg_pwrap(const struct vec3 *r1, const struct vec3 *r2, const struct domain* domain_section);

int is_in_volume(const struct vec3 *vec, const struct domain* domain_section);

#endif //VEC3_MATHS_H
