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
    Recursively fills the octrees of the whole tree with nodes, till there is only one (or zero) node per bottom octree.
    @param tree: Octree to add node into.
    @param tangle Tangle that holds vortices (nodes) to add into the tree (layer of the octree).
*/
void octree_inner_sort(struct octree* tree, const struct tangle_state* tangle)
{
    if (!tree || tree->N <= 1) return; // There is only one node, we do not need to split further.

    octree_make_child_boxes(tree); // Make 8 child boxes.

    for (int k = 0; k < tree->N; ++k) // Loop through all the nodes (vortices) in the current tree.
    {
        for (octree_child_idx child = 0; child < OCTREE_CHILDREN_N; child++)
        {
            // Here we are checking if a given node is in a given child box, continue if yes.
            if (in_octree_box(&tree->children[child]->box, &tangle->vnodes[tree->node_ids[k]]))
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
    tree->box = tangle->box;

    // Sort the points in a recursive function, update physics.
    octree_inner_sort(tree, tangle);
    octree_update_means(tree, tangle);
    return tree;
}

/**
    Checks if the point given by vec is inside the box.
    @param box: Box given by two vectors.
    @param vec: Position of point to check.
    @return Returns True if vec is inside the box, otherwise returns False.
*/
int in_octree_box(const struct domain_box *box, const struct vec3 *vec)
{
    double vx = vec->p[0];
    double vy = vec->p[1];
    double vz = vec->p[2];

    int in_box = 0;

    // Checks if the given point is inside the given box.
    in_box =
        vx >= box->bottom_left_back.p[0] && vx <= box->top_right_front.p[0] &&
        vy >= box->bottom_left_back.p[1] && vy <= box->top_right_front.p[1] &&
        vz >= box->bottom_left_back.p[2] && vz <= box->top_right_front.p[2];

    return in_box;
}

/**
    Helper for creating a child box. Gets the base dimensions from parent_box and returns the box
    corresponding to the child_idx.
    @param parent_box: Parent box is given by two vectors.
    @param child_idx: Id of child box to calculate base dimensions.
    @return Returns base dimensions of child box corresponding to child_idx.
 */
struct domain_box make_child_box(const struct domain_box* parent_box, octree_child_idx child_idx)
{
    //box center and lengths
    struct vec3 c;
    struct vec3 L;

    vec3_add(&c, &parent_box->bottom_left_back, &parent_box->top_right_front);
    vec3_mul(&c, &c, 0.5);

    vec3_sub(&L, &parent_box->top_right_front, &parent_box->bottom_left_back);
    struct domain_box child_box;

    switch (child_idx)
    {
        //back->front -- increasing x
        //left->right -- increasing y
        //bottom->top -- increasing z
    case BLF: //bottom-left-front
        child_box.bottom_left_back = parent_box->bottom_left_back;
        child_box.bottom_left_back.p[0] += L.p[0] / 2;
        child_box.top_right_front = c;
        child_box.top_right_front.p[0] += L.p[0] / 2;
        break;
    case BLB: //bottom-left-back
        child_box.bottom_left_back = parent_box->bottom_left_back;
        child_box.top_right_front = c;
        break;
    case BRF: //bottom-right-front
        child_box.bottom_left_back = c;
        child_box.bottom_left_back.p[2] -= L.p[2] / 2;
        child_box.top_right_front = parent_box->top_right_front;
        child_box.top_right_front.p[2] -= L.p[2] / 2;
        break;
    case BRB: //bottom-right-back
        child_box.bottom_left_back = parent_box->bottom_left_back;
        child_box.bottom_left_back.p[1] += L.p[1] / 2;
        child_box.top_right_front = c;
        child_box.top_right_front.p[1] += L.p[1] / 2;
        break;
    case TLF: //top-left-front
        child_box.bottom_left_back = c;
        child_box.bottom_left_back.p[1] -= L.p[1] / 2;
        child_box.top_right_front = parent_box->top_right_front;
        child_box.top_right_front.p[1] -= L.p[1] / 2;
        break;
    case TLB: //top-left-back
        child_box.bottom_left_back = parent_box->bottom_left_back;
        child_box.bottom_left_back.p[2] += L.p[2] / 2;
        child_box.top_right_front = c;
        child_box.top_right_front.p[2] += L.p[2] / 2;
        break;
    case TRF: //top-right-front
        child_box.bottom_left_back = c;
        child_box.top_right_front = parent_box->top_right_front;
        break;
    case TRB: //top-right-back
        child_box.bottom_left_back = c;
        child_box.bottom_left_back.p[0] -= L.p[0] / 2;
        child_box.top_right_front = parent_box->top_right_front;
        child_box.top_right_front.p[0] -= L.p[0] / 2;
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
        struct domain_box child_box = make_child_box(&tree->box, child);
        tree->children[child] = octree_init(Ninit);
        tree->children[child]->box = child_box;
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
            struct segment seg = seg_pwrap(&tangle->vnodes[idx], &tangle->vnodes[forward], &tangle->box);
            struct vec3 ds = segment_to_vec(&seg);
            vec3_add(&tree->total_circulation, &tree->total_circulation, &ds);
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
    double Lx = tree->box.top_right_front.p[0] - tree->box.bottom_left_back.p[0];
    double Ly = tree->box.top_right_front.p[1] - tree->box.bottom_left_back.p[1];
    double Lz = tree->box.top_right_front.p[2] - tree->box.bottom_left_back.p[2];
    double Lm = Lx > Ly ? (Lx > Lz ? Lx : Lz) : (Ly > Lz ? Ly : Lz); //maximum size

    // Find the distance from r to the box's center of mass.
    struct segment seg = seg_pwrap(r, &tree->centre_of_mass, &tree->tangle->box);
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
