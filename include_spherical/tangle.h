#ifndef TANGLE_H
#define TANGLE_H

#include <stdlib.h>
#include <stdint.h>

#include "vec3_math.h"
#include "boundaries.h"
#include "octree.h"

#define PERIODIC_X 1
#define PERIODIC_Y 1<<1
#define PERIDOIC_Z 1<<2

/*******************************************************************************
 ******************************** TANGLE ***************************************
 *******************************************************************************/

 /*
	The structure that holds the id of next and previous points..
 */
struct neighbour_t {
	int forward;
	int reverse;
};

/*
	The structure that holds all the tangle information.
 */
struct tangle_state {
	struct vec3 *vnodes;			// Positions of the nodes.
	struct vec3 *vels;				// Node velocities. From this we calculate next positions.
	struct vec3 *vs;				// Superfluid velocity at the node.
	struct vec3 *tangents;			// Tangents to the vortices at nodes.
	struct vec3 *normals;			// Normals to the vortices at nodes.

	int *recalculate;				// Flags that the properties of the nodes, currently used for not reconnecting twice in a single pass.

	struct neighbour_t *connections;// Remembers the next and previous points.

	struct node_status *status;		// Status of the node: EMPTY, FREE, PINNED, PINNED_SLIP.

	struct domain domain_section;			// The size of the domain box and the boundary conditions in the 6 cardinal directions (OPEN, PERIODIC, MIRROR).

	int N;							// Total number of all nodes in the tangle.
	int next_empty;					// Id of next empty node in the tangle. We use it to not always search for empty points. Empty point is usually the next one.
	int total_empty;				// Total number of empty nodes in the tangle.
};

/*******************************************************************************
 *************************** TANGLE BASIC FUNCTIONS ****************************
 ******************************************************************************/

void create_tangle(struct tangle_state *tangle, size_t n);

void expand_tangle(struct tangle_state *tangle, size_t n);

void free_tangle(struct tangle_state *tangle);

/*******************************************************************************
 ************************** TANGLE POINTS FUNCTIONS ****************************
 ******************************************************************************/

int number_of_empty_nodes(struct tangle_state *tangle);

int number_of_used_nodes(struct tangle_state *tangle);

struct vec3 step_node(const struct tangle_state *tangle, int i, int where);

int get_next_empty_node(struct tangle_state *tangle);
//add a point between p and p+1 (p+1 in the sense of connections)
int add_point(struct tangle_state *tangle, int p);

void remove_point(struct tangle_state *tangle, int point_idx);

/*******************************************************************************
 ********************** TANGLE FUNCTIONs FOR SIMULATION ************************
 ******************************************************************************/

void enforce_periodic_boundaries(struct tangle_state* tangle);

void enforce_wall_boundaries(struct tangle_state* tangle);

void update_tangent_normal(struct tangle_state *tangle, size_t k);

void update_tangents_normals(struct tangle_state *tangle);

void update_velocity(struct tangle_state* tangle, struct tangle_state* inner_img, struct tangle_state* outer_img, int k, double t,
	struct octree* tree, struct octree* inner_img_tree, struct octree* outer_img_tree);

void update_one_velocity(struct tangle_state* tangle,  int k, double t);

void update_velocities(struct tangle_state *tangle, double t);

void remesh(struct tangle_state *tangle, double min_dist, double max_dist);

void eliminate_small_loops(struct tangle_state *tangle);

static inline void update_tangle(struct tangle_state *tangle, double t)
{
  update_tangents_normals(tangle);
  update_velocities(tangle, t);
}

#endif //TANGLE_H
