#include <math.h>
#include <stdio.h>

#include "vec3_math.h"
#include "boundaries.h"
#include "vortex_constants.h"


/*******************************************************************************
 ********************** PERIODIC AND MIRROR GEOMETRIES *************************
 ******************************************************************************/

 /*
	Inward-facing normals of the domain surface walls.
 */
struct vec3 surface_normal(const struct vec3* where, int surface)
{
	struct vec3 out = *where;
	vec3_normalize(&out);
	
	if (surface == INNER) return out;

	vec3_mul(&out, &out, -1);
	if (surface == OUTER) return out;

	// Suppres warnings.
	return out;
}

 /*
	Periodically wraps the point r2 to follows after the point r1 inside the box.
	In the other words r2 is moved so, that it is in the same period of box as r1.
	@param r1: The reference point.
	@param r2: The next point we want to wrap periodically on the right place.
	@param box: Domain section determining the one period.
	@returns Returns segmen from r1 to r2 that is smaller then one perod of domain box.
 */
struct segment seg_pwrap(const struct vec3 *r1, const struct vec3 *r2, const struct domain *domain_section)
{
	// Angles of the domain.
	double domain_azimut_angle = domain_section->azimut_angle;
    double domain_polar_angle = domain_section->polar_angle;

	// Don't forget that it may be needed to wrap one point multiple times periodically.

	// Check azimutal angles.
	double azimut1 = azimut_angle(r1);
	double azimut2 = azimut_angle(r2);
	while (azimut2 - azimut1 > domain_azimut_angle) azimut2 -= 2 * domain_azimut_angle;
	while (azimut2 - azimut1 < -domain_azimut_angle) azimut2 += 2 * domain_azimut_angle;

	// Check polar angles.
	double polar1 = polar_angle(r1);
	double polar2 = polar_angle(r2);
	while (polar2 - polar1 > domain_polar_angle) polar2 -= 2 * domain_polar_angle;
	while (polar2 - polar1 < -domain_polar_angle) polar2 += 2 * domain_polar_angle;

	// Reverse coordinates.
	double rad2 = radius(r2);

	// And return it as segment.
	struct segment seg = {
		.r1 = *r1,
		.r2 = spherical_to_vector(rad2, azimut2, polar2)
	};

    return seg;
}

/*
	Returns 1 if the given point vec is inside the domain section, otherwise returns 0.
	This function is used in add_point, to see if added point is inside the domain section or not.
	@param vec: Vector to check.
	@param box: Domain section determining the one period.
	@returns Returns 1 if the given point vec is inside the domain section, otherwise returns 0.
*/
int is_in_volume(const struct vec3* vec, const struct domain* domain_section)
{
	// Angles of the domain.
	double domain_azimut_angle = domain_section->azimut_angle;
	double domain_polar_angle = domain_section->polar_angle;

	// Angles of the vector.
	double azimut = azimut_angle(vec);
	double polar = azimut_angle(vec);

	// Returns 1 if inside
	int inside = 1;

	if (fabs(azimut) > domain_azimut_angle || fabs(polar) > domain_polar_angle) inside = 0;

	return inside;
}
