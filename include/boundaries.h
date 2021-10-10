#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "vec3_math.h"

/*******************************************************************************
 ******************************** DOMAIN BOX ***********************************
 ******************************************************************************/

 /*
	Inward-facing normals of the box boundary face walls.
 */
extern const struct vec3 boundary_normals[6];

/*
	Eumerator of wall types (WALL_OPNE, WALL_PERIODIC, WALL_MIRROR).
*/
typedef enum _wt {
    WALL_OPEN, // open wall
    WALL_PERIODIC, // periodic wall
    WALL_MIRROR // mirror wall
} wall_type;

/*
	Enumberator of boundary faces: NOT_A_FACE (-1), LEFT (-x), RIGHT (+x), BACK (-y), FRONT (+y), DOWN (-z), UP (+z).
*/
typedef enum _fc {
    NOT_A_FACE = -1, // not a face
    LEFT, // -x
	RIGHT, // +x
    BACK, // -y
	FRONT, // +y
    DOWN, //-z
	UP // +z
} boundary_faces;

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
  boundary_faces pin_wall; // On which wall the node is pinned: NOT_A_FACE (-1), LEFT (-x), RIGHT (+x), BACK (-y), FRONT (+y), DOWN (-z), UP (+z).
};

/*
	Size of the domain box (bottom_left_back, top_right_front) and its wall types (WALL_OPNE, WALL_PERIODIC, WALL_MIRROR).
*/
struct domain_box {
    struct vec3 bottom_left_back;
    struct vec3 top_right_front;
    wall_type wall[6];
};

/**
	Structure holding information where the image tangle is located and if it should be refleted.
*/
struct image_tangle {
  int shift[3]; // positive/negative number of shifts in the units of box size
  int reflect; // either a boundary_faces enum index or -1 for periodic
};

/**
	Structure holding information about all the image tangles.
	Velocity is then calculated by moving and mirroring the point where we want to find velocity using the information in the image tangle.
	Therefore we don't need to create new tangle containing all the points, we just move the position where we want to calculate the velocity.
*/
struct boundary_images {
  const struct image_tangle *images;
  int n;
};

int in_box(const struct domain_box *box, const struct vec3 *vec);

/*******************************************************************************
 ********************** PERIODIC AND MIRROR GEOMETRIES *************************
 ******************************************************************************/

struct segment seg_pwrap(const struct vec3 *r1, const struct vec3 *r2, const struct domain_box *box);

struct vec3 box_shift(const struct vec3 *v, const struct domain_box *box, int shift[3]);

struct vec3 mirror_shift(const struct vec3 *v, const struct domain_box *box, boundary_faces wall);

struct vec3 mirror_vector_reflect(const struct vec3 *v, boundary_faces wall);

struct vec3 periodic_shift(const struct vec3 *v, const struct domain_box *box, boundary_faces wall);

#endif //VEC3_MATHS_H
