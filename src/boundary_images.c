/**
	Code for initialization and injecion of vorticies.
	@file vortex_injection.c
*/

#include "tangle.h"

/*
	Constant boundary configuration declared in tangle.h.
 */

 // Image tangle configuration of 6 periodic boxes (every wall).
const struct image_tangle periodic_6_img[] = {
    {{-1, 0, 0}, -1},
    {{0, -1, 0}, -1},
    {{0, 0, -1}, -1},
    {{0, 0, 1}, -1},
    {{0, 1, 0}, -1},
    {{1, 0, 0}, -1}
};

// Image tangle configuration of 18 periodic boxes (every wall and edge).
const struct image_tangle periodic_18_img[] = {
    {{-1, -1, 0}, -1},
    {{-1, 0, -1}, -1},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, 0}, -1},
    {{0, -1, -1}, -1},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, -1},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, -1},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, 0}, -1},
    {{1, 0, -1}, -1},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, 0}, -1}
};

// Image tangle configuration of 26 periodic boxes (every wall, edge and corner).
const struct image_tangle periodic_26_img[] = {
    {{-1, -1, -1}, -1},
    {{-1, -1, 0}, -1},
    {{-1, -1, 1}, -1},
    {{-1, 0, -1}, -1},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, -1}, -1},
    {{-1, 1, 0}, -1},
    {{-1, 1, 1}, -1},
    {{0, -1, -1}, -1},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, -1},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, -1},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, -1}, -1},
    {{1, -1, 0}, -1},
    {{1, -1, 1}, -1},
    {{1, 0, -1}, -1},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, -1}, -1},
    {{1, 1, 0}, -1},
    {{1, 1, 1}, -1}
};

// Image tangle configuration periodic in x, y (RIGHT, LEFT, FRONT, BACK) and with walls in z (UP, DOWN).
const struct image_tangle wall_2_4_img[] = {
    {{-1, 0, 0}, -1}, //periodic in x
    {{1, 0, 0}, -1},
    {{0, -1, 0}, -1}, //periodic in y
    {{0, 1, 0}, -1},
    {{0, 0, -1}, DOWN}, //lower z-wall
    {{0, 0, 1}, UP} //upper z-wall
};

// Image tangle configuration periodic in x (RIGHT, LEFT) and with walls in y, z (FRONT, BACK, UP, DOWN).
const struct image_tangle wall_2_2_img[] = {
    {{-1, 0, 0}, -1}, //periodic in x
    {{1, 0, 0}, -1},
    {{0, -1, 0}, BACK}, //front y wall
    {{0, 1, 0}, FRONT},  //back y wall
    {{0, 0, -1}, DOWN}, //lower z-wall
    {{0, 0, 1}, UP} //upper z-wall
};

// Image tangle configuration of parallel planes (walls in DUWN and UP directions).
const struct image_tangle wall_parallel_planes_img[] = {
	{{-1, 0, 0}, -1}, //periodic in x
	{{1, 0, 0}, -1},
	{{0, -1, 0}, -1}, //periodic in y
	{{0, 1, 0}, -1},
	{{0, 0, -1}, DOWN}, //lower z-wall
	{{0, 0, 1}, UP}, //upper z-wall

	{{-1, -1, 0}, -1},
	{{1, -1, 0}, -1},
	{{-1, 1, 0}, -1},
	{{1, 1, 0}, -1},

	{{1, 0, -1}, DOWN},
	{{-1, 0, -1}, DOWN},
	{{0, 1, -1}, DOWN},
	{{0, -1, -1}, DOWN},

	{{1, 0, 1}, UP},
	{{-1, 0, 1}, UP},
	{{0, 1, 1}, UP},
	{{0, -1, 1}, UP},

	{{1, 1, -1}, DOWN},
	{{-1, 1, -1}, DOWN},
	{{-1, 1, -1}, DOWN},
	{{-1, -1, -1}, DOWN},

	{{1, 1, 1}, UP},
	{{-1, 1, 1}, UP},
	{{-1, 1, 1}, UP},
	{{-1, -1, 1}, UP}
};

// Image tangle configuration of canal (dirction of canal is RIGHT to LEFT).
const struct image_tangle wall_canal_img[] = {
	{{-1, 0, 0}, -1}, //periodic in x
	{{1, 0, 0}, -1},
	{{0, -1, 0}, BACK}, //front y wall
	{{0, 1, 0}, FRONT},  //back y wall
	{{0, 0, -1}, DOWN}, //lower z-wall
	{{0, 0, 1}, UP}, //upper z-wall

	{{-1, -1, 0}, BACK},
	{{1, -1, 0}, BACK},
	{{-1, 1, 0}, FRONT},
	{{1, 1, 0}, FRONT},
	{{1, 0, -1}, DOWN},
	{{-1, 0, -1}, DOWN},
	{{1, 0, 1}, UP},
	{{-1, 0, 1}, UP},
};

// Image tangle configuration open everywhere, exept with wall DOWN.
const struct image_tangle wall_1_open_img[] = {
    {{0, 0, -1}, DOWN}
};

// Image tangle configuration periodic in RIGHT, LEFT, FRONT, BACK, UP and with wall DOWN.
const struct image_tangle wall_1_6_img[] = {
    {{-1, 0, 0}, -1},
    {{0, -1, 0}, -1},
    {{0, 0, -1}, DOWN},
    {{0, 1, 0}, -1},
    {{1, 0, 0}, -1}
};

// Image tangle configuration periodic in RIGHT, LEFT, FRONT, BACK, UP plus all edges, and with wall DOWN.
const struct image_tangle wall_1_18_img[] = {
    {{-1, -1, 0}, -1},
    {{-1, 0, -1}, DOWN},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, 0}, -1},
    {{0, -1, -1}, DOWN},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, DOWN},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, DOWN},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, 0}, -1},
    {{1, 0, -1}, DOWN},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, 0}, -1},
};

// Image tangle configuration periodic in RIGHT, LEFT, FRONT, BACK, UP plus all edges and corners, and with wall DOWN.
const struct image_tangle wall_1_26_img[] = {
    {{-1, -1, -1}, DOWN},
    {{-1, -1, 0}, -1},
    {{-1, -1, 1}, -1},
    {{-1, 0, -1}, DOWN},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, -1}, DOWN},
    {{-1, 1, 0}, -1},
    {{-1, 1, 1}, -1},
    {{0, -1, -1}, DOWN},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, DOWN},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, DOWN},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, -1}, DOWN},
    {{1, -1, 0}, -1},
    {{1, -1, 1}, -1},
    {{1, 0, -1}, DOWN},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, -1}, DOWN},
    {{1, 1, 0}, -1},
    {{1, 1, 1}, -1}
};

// Constant boundary configuration open everywhere.
const struct boundary_images open_boundaries = {
    .images = NULL,
    .n = 0
};

// Constant boundary configuration of 6 periodic boxes (every wall).
const struct boundary_images periodic_6 = {
   .images = periodic_6_img,
   .n = 6
};

// Constant boundary configuration of 18 periodic boxes (every wall and edge).
const struct boundary_images periodic_18 = {
   .images = periodic_18_img,
   .n = 18
};

// Constant boundary configuration of 26 periodic boxes (every wall, edge and corner).
const struct boundary_images periodic_26 = {
   .images = periodic_26_img,
   .n = 26
};

// Constant boundary configuration open everywhere, exept with wall DOWN.
const struct boundary_images wall_1_open = {
   .images = wall_1_open_img,
   .n = 1
};

// Constant boundary configuration periodic in RIGHT, LEFT, FRONT, BACK, UP and with wall DOWN.
const struct boundary_images wall_1_6 = {
   .images = wall_1_6_img,
   .n = 6
};

// Constant boundary configuration periodic in RIGHT, LEFT, FRONT, BACK, UP plus all edges, and with wall DOWN.
const struct boundary_images wall_1_18 = {
   .images = wall_1_18_img,
   .n = 18
};

// Constant boundary configuration periodic in RIGHT, LEFT, FRONT, BACK, UP plus all edges and corners, and with wall DOWN.
const struct boundary_images wall_1_26 = {
   .images = wall_1_26_img,
   .n = 26
};

// Constant boundary configuration periodic in x, y (RIGHT, LEFT, FRONT, BACK) and with walls in z (UP, DOWN).
const struct boundary_images wall_2_4 = {
    .images = wall_2_4_img,
    .n = 6
};
 
// Constant boundary configuration periodic in x (RIGHT, LEFT) and with walls in y, z (FRONT, BACK, UP, DOWN).
const struct boundary_images wall_2_2 = {
    .images = wall_2_2_img,
    .n = 6
};

// Constant boundary configuration of parallel planes (walls in DUWN and UP directions).
const struct boundary_images wall_parallel_planes = {
	.images = wall_parallel_planes_img,
	.n = 26
};

// Constant boundary configuration of canal (dirction of canal is RIGHT to LEFT).
const struct boundary_images wall_canal = {
	.images = wall_canal_img,
	.n = 14
};
