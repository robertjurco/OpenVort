#include <math.h>
#include <string.h>
#include <stdio.h>

#include "vec3_math.h"

#include "external_velocity.h"

/*
	Calculates the NOFLOW velocity and saves it into *res.
	@param res: Pointer to vector where the result is stored.
 */
void get_v_noflow(struct vec3 *res)
{
    *res = vec3(0, 0, 0); // No flow.
}

/*
	Calculates the SIMPLE velocity with given strength and direction, and saves it into *res.
	@param res: Pointer to vector where the result is stored.
	@param dir: Direction of SIMPLE flow.
	@param strength: Strength of SIMPLE flow [cm/s].
 */
void get_v_simple(struct vec3 *res, struct vec3 *dir, double strength)
{
    *res = *dir;
    vec3_normalize(res);
    vec3_mul(res, res, strength);
}

/*
	Calculates the DECREASING velocity with given starting strength, final strength, time of linear decrease and direction, and saves it into *res.
	@param res: Pointer to vector where the result is stored.
	@param dir: Direction of DECREASING flow.
	@param starting_strength: Starting strength for the DECREASING flow type.
	@param time_of_linear_decrease: Time for which the DECREASING flow type is lineary decreasing from its starting strength onto strength.
	@param strength: Final constant strength for the DECREASING flow type.
 */
void get_v_decreasing(struct vec3 *res, double t,struct vec3 *dir, double starting_strength, double time_of_liear_decrease, double strength)
{
	double str;

    if (t < time_of_liear_decrease) {
        str = starting_strength - t * (starting_strength - strength) / time_of_liear_decrease;
    } else {
		str = strength;
    }

    *res = *dir;
    vec3_normalize(res);
    vec3_mul(res, res, str);
}

/*
	Calculates the fluid velocity at point *where and saves it into *res. This function is called by public get_vn, get_vb, get_vs functions.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@param ext_vel_param: Parameter of external velocity. See ext_vel_param in file external_velocity.h.
	@returns Returns 1 on success.
 */
int get_v(const struct vec3 *where, double t, struct vec3 *res, struct ext_vel_param *ext_vel_param)
{
    switch (ext_vel_param->type) {
        case NO_FLOW:
            get_v_noflow(res);
            return 1;
        case SIMPLE:
            get_v_simple(res, &ext_vel_param->direction, ext_vel_param->strength);
            return 1;
        case DECREASING:
            get_v_decreasing(res, t, &ext_vel_param->direction, ext_vel_param->starting_strength, ext_vel_param->time_of_linear_decrease, ext_vel_param->strength);
            return 1;
        default:
            printf("Unknown external velocity type.");
            return 0;
    }
}

/*
	Calculates the normal fluid velocity at point *where and time t, and saves it into *res.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@returns Returns 1 on success.
 */
int get_vn(const struct vec3 *where, double t, struct vec3 *res)
{
  return get_v(where, t, res, &vn_conf);
}

/*
	Calculates the superfluid velocity at point *where and time t, and saves it into *res.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@returns Returns 1 on success.
 */
int get_vs(const struct vec3 *where, double t, struct vec3 *res)
{
  return get_v(where, t, res, &vs_conf);
}

/*
	Calculates the boundary velocity at point *where and time t, and saves it into *res.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@returns Returns 1 on success.
 */
int get_vb(const struct vec3 *where, double t, struct vec3 *res)
{
  return get_v(where, t, res, &vb_conf);
}
