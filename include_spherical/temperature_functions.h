#ifndef TEMPERATURE_FUNCTIONS_H
#define TEMPERATURE_FUNCTIONS_H

#include "vec3_math.h"

extern double temp_la;

double total_density(double T);

double normal_density(double T);

double superfluid_density(double T);

double entropy(double T);

double mutual_firction_alpha(double T);

double mutual_firction_alpha_p(double T);

/*
	Calculates the normal fluid velocity at point *where and time t, and saves it into *res.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@returns Returns 1 on success.
 */
struct vec3 get_vn(const struct vec3* where, double t);

/*
	Calculates the superfluid velocity at point *where and time t, and saves it into *res.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@returns Returns 1 on success.
 */
struct vec3 get_vs(const struct vec3* where, double t);

double get_temperature(const struct vec3* where);

#endif //EMPERATURE_FUNCTIONS_H
