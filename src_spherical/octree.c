/**
    Code for Banes-Hut tree algorithm.
    @file octree.c

    Implementation:

    struct octree *tree = NULL;
    if(use_BH) tree = octree_build(tangle);
    octree_get_vs(tangle, tree, &tangle->vnodes[k], BH_resolution, &velocity);
    octree_destroy(tree);
*/

#include <stdlib.h>
#include <math.h>

#include "stdio.h"
#include "vortex_constants.h"
#include "octree.h"
#include "vec3_math.h"

/**
    Initialize octree (layer of the tree) and array of nodes (vortices) inside the octree.
    @param Ninit: Number of nodes (vortices) inside the current octree.
    @return The initialized octree.
*/
struct octree* octree_init(int Ninit)
{
    struct octree *tree = (struct octree*)malloc(sizeof(struct octree));

    tree->N=0;
    tree->N_total = Ninit;
    tree->node_ids = (int*)malloc(sizeof(int)*Ninit);

    return tree;
}

/**
    Destroyes (frees) the tree structure.
    @param tree: Octree to free.
*/
void octree_destroy(struct octree *tree)
{
    if(!tree) return;

    //we need to destroy it from the bottom up
    for(octree_child_idx child = 0; child < OCTREE_CHILDREN_N; child++)
    {
         if(tree->children[child]) octree_destroy(tree->children[child]);
    }
    free(tree->node_ids);
    free(tree);
}

/**
    Adds a node to the root of the octree, this DOES NOT add the node recursively.
    @param tree: Octree to add node in.
    @param node_id: ID od node (vortice) to add into the tree (layer of the octree).
    @return Returns current octree if node is sucesfully added, NULL if not.
*/
struct octree *octree_add_node(struct octree *tree, int node_id)
{
    //adds a node to the root of the tree
    //this DOES NOT add the node recursively
    if(tree->N < tree->N_total)
    {
        // increment N
        tree->node_ids[tree->N++] = node_id;
        return tree;
    }

    tree->node_ids = (int*)realloc(tree->node_ids, sizeof(int)*2*tree->N_total);
    if(tree->node_ids)
    {
        tree->N_total *= 2;
        // increment N
        tree->node_ids[tree->N++] = node_id;
        return tree;
    }
    else return NULL;
}

/**
	Checks if the point given by vec is inside the box.
	@param box: Box given by two vectors.
	@param vec: Position of point to check.
	@return Returns True if vec is inside the box, otherwise returns False.
*/
int in_octree_box(const struct tree_domain* domain_section, const struct vec3* vec)
{
	double r = radius(vec);
	double az = azimut_angle(vec);
	double pl = polar_angle(vec);

	int in_box = 0;

	// Checks if the given point is inside the given box.
	in_box =
		r  >= domain_section->inner_radius	  && r  <= domain_section->outer_radius    &&
		az >= domain_section->negative_azimut && az <= domain_section->positive_azimut &&
		pl >= domain_section->negative_polar  && pl <= domain_section->positive_polar;

	return in_box;
}

/**
    Recursively fills the octrees of the whole tree with nodes, till there is only one (or zero) node per bottom octree.
    @param tree: Octree to add node into.
    @param tangle Tangle that holds vortices (nodes) to add into the tree (layer of the octree).
*/
void octree_inner_sort(struct octree* tree, const struct tangle_state* tangle)
{
    if (!tree || tree->N <= 1) return; // There is only one node, we do not need to split further.

	/*
	printf("\nri %g, ro %g, pa %g, na %g, pp %g, np %g\n", tree->domain_section.inner_radius, tree->domain_section.outer_radius, 180 / M_PI *tree->domain_section.positive_azimut, 180 / M_PI * tree->domain_section.negative_azimut, 180 / M_PI * tree->domain_section.positive_polar, 180 / M_PI * tree->domain_section.negative_polar);
	fflush(stdout);
	double rrr = radius(&tangle->vnodes[tree->node_ids[3]]);
	double aaa = 180 / M_PI * azimut_angle(&tangle->vnodes[tree->node_ids[1]]);
	double ppp = 180 / M_PI * polar_angle(&tangle->vnodes[tree->node_ids[1]]);
	printf("r %g, a %g, p %g\n", rrr, aaa, ppp);
	fflush(stdout);
	*/

    octree_make_child_boxes(tree); // Make 8 child boxes.

    for (int k = 0; k < tree->N; ++k) // Loop through all the nodes (vortices) in the current tree.
    {
        for (octree_child_idx child = 0; child < OCTREE_CHILDREN_N; child++)
        {
            // Here we are checking if a given node is in a given child box, continue if yes.
            if (in_octree_box(&tree->children[child]->domain_section, &tangle->vnodes[tree->node_ids[k]]))
            {
                // Add the node into right child box, then break the loop.
                if (!octree_add_node(tree->children[child], tree->node_ids[k])) printf("failed to add node");
                break;
            }
        }
    }

    // Continue recursively till the bottom of the tree.
    for (octree_child_idx child = 0; child < OCTREE_CHILDREN_N; child++) octree_inner_sort(tree->children[child], tangle);
}

/**
    Builds the whole tree from vortices (nodes) of tangle.
    @param tangle: Tangle that holds vortices (nodes) to make the tree from.
    @return Returns top octree of the builded tree.
*/
struct octree *octree_build(const struct tangle_state *tangle)
{
    struct octree *tree = octree_init(tangle->N);
    //load up all the non-empty points

    for(int k=0; k < tangle->N; ++k)
    {
        if(tangle->status[k].status == EMPTY) continue;
        if(!octree_add_node(tree, k)) printf("ERROR: Failed to add node while building the octree.\n");
    }

    tree->tangle = tangle;

	tree->domain_section.outer_radius    =  tangle->domain_section.outer_radius;
	tree->domain_section.inner_radius    =  tangle->domain_section.inner_radius;
	tree->domain_section.positive_azimut =  tangle->domain_section.azimut_angle;
	tree->domain_section.negative_azimut = -tangle->domain_section.azimut_angle;
	tree->domain_section.positive_polar  =  tangle->domain_section.polar_angle;
	tree->domain_section.negative_polar  = -tangle->domain_section.polar_angle;

    // Sort the points in a recursive function, update physics.
    octree_inner_sort(tree, tangle);

    octree_update_means(tree, tangle);

    return tree;
}

/**
    Helper for creating a child box. Gets the base dimensions from parent_box and returns the box
    corresponding to the child_idx.
    @param child_idx: Id of child box to calculate base dimensions.
    @return Returns base dimensions of child box corresponding to child_idx.
 */
struct tree_domain make_child_box(const struct tree_domain* parent_section, octree_child_idx child_idx)
{
	double inner  = parent_section->inner_radius;
	double outer  = parent_section->outer_radius;
	double pos_az = parent_section->positive_azimut;
	double neg_az = parent_section->negative_azimut;
	double pos_pl = parent_section->positive_polar;
	double neg_pl = parent_section->negative_polar;

    struct tree_domain child_box;

    switch (child_idx)
    {
        //back->front -- increasing x
        //left->right -- increasing y
        //bottom->top -- increasing z
    case BLF: //bottom-left-front
		child_box.inner_radius    = (outer + inner) / 2;
		child_box.outer_radius    = outer;
		child_box.positive_azimut = (pos_az + neg_az) / 2;
		child_box.negative_azimut = neg_az;
		child_box.positive_polar  = (pos_pl + neg_pl) / 2;
		child_box.negative_polar  = neg_pl;
        break;
    case BLB: //bottom-left-back
		child_box.inner_radius    = inner;
		child_box.outer_radius    = (outer + inner) / 2;
		child_box.positive_azimut = (pos_az + neg_az) / 2;
		child_box.negative_azimut = neg_az;
		child_box.positive_polar  = (pos_pl + neg_pl) / 2;
		child_box.negative_polar  = neg_pl;
        break;
    case BRF: //bottom-right-front
		child_box.inner_radius    = (outer + inner) / 2;
		child_box.outer_radius    = outer;
		child_box.positive_azimut = (pos_az + neg_az) / 2;
		child_box.negative_azimut = neg_az;
		child_box.positive_polar  = pos_pl;
		child_box.negative_polar  = (pos_pl + neg_pl) / 2;
        break;
    case BRB: //bottom-right-back
		child_box.inner_radius    = inner;
		child_box.outer_radius    = (outer + inner) / 2;
		child_box.positive_azimut = (pos_az + neg_az) / 2;
		child_box.negative_azimut = neg_az;
		child_box.positive_polar  = pos_pl;
		child_box.negative_polar  = (pos_pl + neg_pl) / 2;
        break;
    case TLF: //top-left-front
		child_box.inner_radius    = (outer + inner) / 2;
		child_box.outer_radius    = outer;
		child_box.positive_azimut = pos_az;
		child_box.negative_azimut = (pos_az + neg_az) / 2;
		child_box.positive_polar  = (pos_pl + neg_pl) / 2;
		child_box.negative_polar  = neg_pl;
        break;
    case TLB: //top-left-back
		child_box.inner_radius    = inner;
		child_box.outer_radius    = (outer + inner) / 2;
		child_box.positive_azimut = pos_az;
		child_box.negative_azimut = (pos_az + neg_az) / 2;
		child_box.positive_polar  = (pos_pl + neg_pl) / 2;
		child_box.negative_polar  = neg_pl;
        break;
    case TRF: //top-right-front
		child_box.inner_radius    = (outer + inner) / 2;
		child_box.outer_radius    = outer; 
		child_box.positive_azimut = pos_az;
		child_box.negative_azimut = (pos_az + neg_az) / 2;
		child_box.positive_polar  = pos_pl;
		child_box.negative_polar  = (pos_pl + neg_pl) / 2;
        break;
    case TRB: //top-right-back
		child_box.inner_radius    = inner;
		child_box.outer_radius    = (outer + inner) / 2;
		child_box.positive_azimut = pos_az;
		child_box.negative_azimut = (pos_az + neg_az) / 2;
		child_box.positive_polar  = pos_pl;
		child_box.negative_polar  = (pos_pl + neg_pl) / 2;
        break;
    default:
        printf("make_child_box: incorrect child index");
        break;
    }

    return child_box;
}

/**
    Creates children and intialize them for the current octree.
    @param tree: Octree to create and initialize children for.
 */
void octree_make_child_boxes(struct octree *tree)
{
    int Ninit = tree->N_total / 8;
    if(Ninit == 0) Ninit = 8; //less than 8 elements

    for(octree_child_idx child = 0; child < OCTREE_CHILDREN_N; ++child)
    {
        struct tree_domain child_box = make_child_box(&tree->domain_section, child);
        tree->children[child] = octree_init(Ninit);
        tree->children[child]->domain_section = child_box;
        tree->children[child]->tangle = tree->tangle;
    }
}

/**
    Updates the center of mass and total circluation for every octree of the whole tree.
    @param tree: Top octree of the tree to update the center of mass and total circluation for.
    @param tangle: Tangle holding informations to calculate the center of mass and total circluation. tree->node_ids[i] ID of i-th point
    @param in the tree (layer of octree), it's position is stored in tangle->vnodes[tree->node_ids[i]]
 */
void octree_update_means(struct octree *tree, const struct tangle_state *tangle)
{
    // Checks if the tree exists (we are not yet on the bottom of the octree).
    if(!tree || tree->N == 0) return;

    // Calculates centre of mass.
    tree->centre_of_mass = vec3(0,0,0);
    for(int k=0; k < tree->N; k++)
    {
        vec3_add(&tree->centre_of_mass, &tree->centre_of_mass, &tangle->vnodes[tree->node_ids[k]]);
    }
    vec3_mul(&tree->centre_of_mass, &tree->centre_of_mass, 1.0/tree->N);

    // Calculates total circulation.
    tree->total_circulation = vec3(0,0,0);
    for(int k=0; k < tree->N; k++)
    {
        int idx = tree->node_ids[k];
        int forward = tangle->connections[idx].forward;
        if (forward > -1) {
            struct segment seg = seg_pwrap(&tangle->vnodes[idx], &tangle->vnodes[forward], &tangle->domain_section);
            struct vec3 ds = segment_to_vec(&seg);
			if (vec3_len(&ds) > 5 * global_dl_max) {
				printf("ERROR: length of segment is > 5 * global_dl_max in octree\n");
			} else {
				vec3_add(&tree->total_circulation, &tree->total_circulation, &ds);
			}
        }
    }

    // Recursively continues till the bottom of the octree.
    for(octree_child_idx child=0; child < OCTREE_CHILDREN_N; child++) octree_update_means(tree->children[child], tangle);
}

/**
    Calculate the velocity induced by whole octree at point r and saves it into vector res.
    @param tangle: Tangle is needed for tree->tangle->box to pwrap vectors from r to points of octree.
    @param tree: Tree structure given by the top octree.
    @param r: Position where we want to calculate induced velocity of superfluid component.
    @param resolution: Resolution of Barnes-Hut approximation: the box width / the distance to the center of mass of box from r.
    @param res: Vector to store output.
 */
void octree_get_vs(const struct tangle_state *tangle, const struct octree *tree, const struct vec3 *r, double resolution, struct vec3 *res)
{
    *res = vec3(0,0,0);

    if(tree->N == 0) return;

    // Find the biggest side of the box.
    double rd = tree->domain_section.outer_radius - tree->domain_section.inner_radius;
	double az = tree->domain_section.outer_radius * 2 * fabs(sin((tree->domain_section.positive_azimut - tree->domain_section.negative_azimut)/2));
    double pl = tree->domain_section.outer_radius * 2 * fabs(sin((tree->domain_section.positive_polar - tree->domain_section.negative_polar)/2)) 
		* fabs(cos((tree->domain_section.positive_azimut + tree->domain_section.negative_azimut) / 2));
    double Lm = rd > az ? (rd > pl ? rd : pl) : (az > pl ? az : pl); //maximum size

    // Find the distance from r to the box's center of mass.
    //struct segment seg = seg_pwrap(r, &tree->centre_of_mass, &tree->tangle->domain_section);
	struct segment seg;
	seg.r1 = *r;
	seg.r2 = tree->centre_of_mass;
    double d = segment_len(&seg);

    // Too close, use LIA.
    if(d < 1e-8) return;

    // Check the resolution, if it's not good enough go deeper to the octree.
    if(Lm/d < resolution || tree->N == 1)
    {
        struct vec3 R; // Vector from r to the box's center of mass.
        struct vec3 Rp1; // Vector from r to the box's (center of mass + total circulation).
        struct vec3 q = tree->total_circulation;
        R = segment_to_vec(&seg);
        vec3_add(&Rp1, &R, &q);

        double lR   = vec3_len(&R); // Distance from r to the box's center of mass.
        double lRp1 = vec3_len(&Rp1); // Distance from r to the box's (center of mass + total circulation).
        double denom = lR*lRp1*(lR*lRp1 + vec3_dot(&R, &Rp1)); // Denominator.
        double f = KAPPA/4/M_PI; // Constant coefficient.

        // This can happen in periodic boundary conditions.
        if(lR < 1e-8 || lRp1 < 1e-8)
          return;

        // If R and Rp1 are colinear, the result is 0, but code below would try to calculate 0/0.
        if(vec3_ndot(&R, &Rp1) < 1e-8)
          return;

        // Calculate velocity.
        vec3_cross(res, &R, &Rp1);
        vec3_mul(res, res, f*(lR + lRp1)/denom);

        return;
    }

    // If resolution is not sufficient, we need to go deeper into the octree.
    struct vec3 partial_vs, total_vs;
    total_vs = vec3(0,0,0);
    for(int kk = 0; kk < 8; ++kk)
    {
        octree_get_vs(tangle, tree->children[kk], r, resolution, &partial_vs);
        vec3_add(&total_vs, &total_vs, &partial_vs);
    }
    *res = total_vs;
}
