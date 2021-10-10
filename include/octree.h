/**
	Header for Banes-Hut tree algorithm.
	@file octree.h

	Implementation:

	struct octree *tree = NULL;
	if(use_BH) tree = octree_build(tangle);
	octree_get_vs(tangle, tree, &tangle->vnodes[k], BH_resolution, &velocity);
	octree_destroy(tree);
*/

#ifndef INCLUDE_OCTREE_H
#define INCLUDE_OCTREE_H

#include "vec3_math.h"
#include "tangle.h"

/**
	Enumerates the child boxes of the octree. NA=-1 represents no children. BLF=0, BLB=1, BRF=2, BRB=3, 
	TLF=4, TLB=5, TRF=6, TRB=7, represents Bottom {Left/Right} {Front/Back} and Top {Left/Right} {Front/Back} children.
	OCTREE_CHILDREN_N=8 represents the total number of children, that is always 8.
*/
typedef enum _octree_child_idx {
  NA=-1,
  BLF, BLB, BRF, BRB, //Bottom {Left/Right} {Front/Back}
  TLF, TLB, TRF, TRB, //Top {Left/Right} {Front/Back}
  OCTREE_CHILDREN_N
} octree_child_idx;


/**
	Structure holding all informations about otree.
	@param N_total: Total number of nodes (representing points of vortices) in tangle.
	@param N: Number of nodes in this octree
	@param node_ids: Ids of nodes (representing points of vortices) in this octree
	@param centre_of_mass: Center of mass of this octree used to calulate induced superfluid velocity.
	@param total_circulation: Total circulation of this octree used to calulate induced superfluid velocity.
	@param box: Box of this octree given by vectors bottom_left_back and top_right_front.
	@param children: Children (sub-octrees) of this octree. They have the same structure as this octree.
	@param tangle: Reference to the tangle used to build the octree.
*/
struct octree {
  int N_total;
  int N;
  int *node_ids;
  struct vec3 centre_of_mass;
  struct vec3 total_circulation;
  struct domain_box box;
  struct octree *children[8];
  const struct tangle_state *tangle;
};

struct octree* octree_init(int Ninit);
void octree_destroy(struct octree *tree);

struct octree *octree_build(const struct tangle_state *tangle);

// helper functions

void octree_make_child_boxes(struct octree *tree);
void octree_update_means(struct octree *tree, const struct tangle_state *tangle);

// velocity calculation

void octree_get_vs(const struct tangle_state *tangle, const struct octree *tree, const struct vec3 *r, double resolution, struct vec3 *res);

#endif /* INCLUDE_OCTREE_H */
