#ifndef VORTEX_UTILS_H
#define VORTEX_UTILS_H

//This file contains miscellaneous utilities related to
//the handling of the tangle. Actual calculations belong to tangle.h.

//This contains routines for creating basic vortex geometries and
//handling file I/O

#include <stdlib.h>
#include <stdio.h>
#include "vec3_math.h"
#include "tangle.h"


/*
 * Debugging utilities
 */


/*
 * Wall-related stuff
 */

 void clip_at_wall(struct tangle_state *tangle);

 double wall_dist(const struct tangle_state *tangle, int k, boundary_faces wall);

// Retutns length of a vortex. Vortex is given by starting point. In case of an open vortes connection forward should be -1.
 double vortex_length(const struct tangle_state *tangle, int starting_point);

#endif
