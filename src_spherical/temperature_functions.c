#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vec3_math.h"
#include "vortex_constants.h"
#include "temperature_functions.h"

double temp_la = 2.1768;

double Spline(int array_length, double* A, double* Q, double T) 
{
	int N = array_length - 5 - 1;

	if (T < Q[4] || T > Q[N + 5]) {
		printf("Temperature out of bounds.");
		exit(EXIT_FAILURE);
		return 0.0;
	}

	// Finds I, J so that Q[I+4] =< T < Q[J+4].
	int I = 0;
	int J = N + 1;
	int K;
	while ((J - I) > 1) {
		K = floor((I + J) / 2);
		if (T >= Q[K + 4]) I = K;
		else  J = K;
	}

	// Copy temperatures from J to J+3.
	int JP3 = J + 3;
	double D[JP3+1];
	for (I = J; I <= JP3; I++) D[I] = A[I];
	
	// Three times recalculate.
	for (K = 1; K <= 3; K++) {
		int JPL = JP3 - K;
		for (I = J; I <= JPL; I++) {
			D[I] = ((T - Q[I + K]) * D[I + 1] + (Q[I + 4] - T) * D[I]) / (Q[I + 4] - Q[I + K]);
		}
	}
	
	double out = D[J];
	return out;
}


double superfluid_density(double T)
{
	printf("entering sup density %g", T);
	if (T >= 0.0 && T <= temp_la) {
		double Nr = 26;
		double Qr[26] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.443, 0.9012, 1.5419, 1.7540, 1.918, 2.111, 2.156991, 2.173218, 2.175647, 2.176358, 2.176568, 2.176692, 2.176766, 2.176791, 2.176798, 2.176799, 2.17679999, 2.1768, 2.1768, 2.1768, 2.1768 };
		double Ar[26] = { 0.0, 1.451275432822459e-1, 1.451334563362309e-1, 1.449759191497576e-1, 1.455008000684433e-1, 1.4075e-1, 1.095e-1, 8.15e-2, 5.30e-2, 2.1e-2, 8.904576e-3, 3.053214e-3, 1.494043e-3, 8.342826e-4, 5.10686e-4, 2.8379e-4, 1.287426e-4, 5.202569e-5, 2.153580e-5, 8.564206e-6, 3.567958e-6, 0.0 };
		double deni = Spline(Nr, Ar, Qr, T);
		printf("finised sup dens");
		return(deni * 1000);
	}
	else {
		printf("Temperature out of bounds in superfluid_density(double T).");
		exit(EXIT_FAILURE);
		return 0.0;
	}
}

double total_density(double T)
{
	if (T >= 0 && T <= 4.9) {
		double ro = 0.0;
		double ro_0 = 0.1451397;
		double ro_la = 0.1461087;
		double delt_ro = 0.0;

		double* a;
		double* b;
		double a1[3] = { 0.0, -7.57537e-3, 6.87483e-3 }; //1.344 - temp_la
		double a2[3] = { 0.0, -7.94605e-3, 5.07051e-3 }; //temp_la - 4.9
		double b1[8] = { 0.0, 3.79937e-3, 1.86557e-3, 4.88345e-3, 0.0, 0.0, 0.0, 0.0 }; //1.344 - temp_la
		double b2[8] = { 0.0, -30.3511e-3, -10.2326e-3, -3.00636e-3, 0.240720e-3, -2.45749e-3, 1.53454e-3, -0.308182e-3 }; //temp_la - 4.9

		double m[5] = { 0.0, -1.26935e-5, 7.12413e-5, -16.7461e-5, 8.75342e-5 };
		double t = T - temp_la;
		ro = 0.0;
		if (T < 1.344) { //eq(1.2)
			for (int i = 1; i <= 4; i++) {
				ro = ro + m[i] * pow(T, i + 1);
			}
			ro = ro + ro_0; //formula for ro
		}
		else if (T >= 1.344 && T <= 4.9) { //eq(1.4)
			if (T <= temp_la) {
				a = a1;
				b = b1;
			}
			else {
				a = a2;
				b = b2;
			}
			delt_ro = 0;
			for (int i = 1; i <= 2; i++) {
				if (t != 0) {
					delt_ro = delt_ro + a[i] * pow(t, i) * log(fabs(t));
				}
				else {
					delt_ro = delt_ro;
				}
			}
			for (int i = 1; i <= 7; i++) {
				delt_ro = delt_ro + b[i] * pow(t, i);
			}
			ro = ro_la * (delt_ro + 1);
		}
		return(ro * 1000); //*1000 for units
	}
	else {
		printf("Temperature out of bounds in total_density(double T).");
		exit(EXIT_FAILURE);
		return 0.0;
	}
}

double normal_density(double T)
{
	if (T >= 1.1 && T <= temp_la)
	{
		double total_den = total_density(T);
		double super_den = superfluid_density(T);
		double deni = total_den - super_den;
		return deni;
	}
	else {
		printf("Temperature out of bounds in normal_density(double T).");
		exit(EXIT_FAILURE);
		return 0.0;
	}
}

double entropy(double T) {
	if (T >= 0.1 && T <= 5.0) {
		int Ns = 30;
		double Qs[30] = { 0.0, 0.10000, 0.10000, 0.10000, 0.10000, 0.11123, 0.21625, 0.27672, 0.44817, 0.56369, 0.71969, 0.85360, 1.00186, 1.40184, 1.72912, 1.99887, 2.07775, 2.15384, 2.17680, 2.17680, 2.17680, 2.21160, 2.34820, 2.69075, 4.38277, 4.76361, 5.00000, 5.0, 5.0, 5.0 };
		double As[30] = { 0.0, 6.60e-6, 7.60e-6, 1.60e-5, 4.61413e-5, 1.8400e-4, 4.68675e-4, 0.00121, 0.00236, 0.00561, 0.01978, 0.08337, 0.35824, 0.76882, 1.17803, 1.38629, 1.53613, 1.58271, 1.62419, 1.72261, 1.91862, 2.53529, 3.16601, 4.10379, 4.44077, 4.61570 };
		double deni = Spline(Ns, As, Qs, T);
		return(deni * 1000);
	}
	else {
		printf("Temperature out of bounds in entropy(double T).");
		exit(EXIT_FAILURE);
		return 0.0;
	}
}

double mutual_firction_alpha(double T)
{
	if (T >= 1.3 && T < temp_la) {
		double denib;

		if (T <= 2.167) {
			double ep = log10(1 - T / temp_la);
			int Ns = 12;
			double Qs[12] = { 0.0, -5.0, -5.0, -5.0, -5.0, -2.5, -2.0, -0.8, -0.387958059947, -0.387958059947, -0.387958059947, -0.387958059947 };
			double As[12] = { 0.0, 1.31928144433, 1.12452707801, 0.639314792565, 0.313383532495, -0.162687403543, 0.092047691284, 0.188452616588 };
			denib = Spline(Ns, As, Qs, ep);
			denib = pow(10, denib);

		}
		else {
			denib = 0.47 * pow(1 - T / temp_la, -0.33);
		}
		double denia = (denib / 2) * (1 - superfluid_density(T) / total_density(T));

		return denia;
	}
	else {
		printf("Temperature out of bounds in mutual_firction_alpha(double T).");
		exit(EXIT_FAILURE);
		return 0.0;
	}
}

double mutual_firction_alpha_p(double T)
{
	if (T >= 1.3 && T < temp_la) {
		double denib2;

		if (T <= 2.134) {
			double ep = log10(1 - T / temp_la);
			double Ns = 13;
			double Qs[13] = { 0.0, -5.0, -5.0, -5.0, -5.0, -3.55, -3.2, -2.5, -1.0, -0.384067377871, -0.384067377871, -0.384067377871, -0.384067377871 };
			double As[13] = { 0.0, -8.47218032526e-2, 0.931621715174, 0.973263359433, 1.10543591819, 1.15904485127, 1.18311634566, 1.17480594214, 1.19458392766 };
			denib2 = Spline(Ns, As, Qs, ep);
			denib2 = pow(10, denib2) - 15;
		}
		else {
			denib2 = -0.34 * pow(1 - T / temp_la, -0.33) + 1.01;
		}
		double denia2 = (denib2 / 2) * (1 - superfluid_density(T) / total_density(T));

		return denia2;
	}
	else {
		printf("Temperature out of bounds in mutual_firction_alpha_p(double T).");
		exit(EXIT_FAILURE);
		return 0.0;
	}
}

/*
	Calculates the normal fluid velocity at point *where and time t, and saves it into *res.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@returns Returns 1 on success.
 */
struct vec3 get_vn(const struct vec3* where, double T)
{
	struct vec3 out = *where;

	vec3_normalize(&out);

	double r = radius(where) * 0.01;
	// v_n = Q_dot / ( T * s * rho * 4 * pi * r^2)
	double v_n = q_dot / (T * entropy(T) * total_density(T) * 4 * M_PI * pow(r, 2));
	
	vec3_mul(&out, &out, v_n * 100);

	return out;
}

/*
	Calculates the superfluid velocity at point *where and time t, and saves it into *res.
	@param where: Where to calculate the velocity.
	@param t: In what time to calculate velocity.
	@param res: Pointer to vector where the result is stored.
	@returns Returns 1 on success.
 */
struct vec3 get_vs(const struct vec3* where, double T)
{
	struct vec3 out = *where;

	vec3_normalize(&out);

	double r = radius(where) * 0.01;

	// v_s = - rho_n * Q_dot / ( T * s * rho_s * rho * 4 * pi * r^2)
	double v_s = - q_dot * normal_density(T) / (T * entropy(T) * superfluid_density(T) * total_density(T) * 4 * M_PI * pow(r, 2));

	vec3_mul(&out, &out, v_s * 100); // To centimeters.

	return out;
}

double gamma(double T)
{
	return (155.8 * T - 75) * 10000; // To meters squared.
	//return (74.9 * T - 70.7) * 10000; // To meters squared.
}

double get_temperature(const struct vec3* where)
{
	double x = radius(where) * 0.01; // To meters.

	double v_n = q_dot / (temp_bath * entropy(temp_bath) * total_density(temp_bath) * 4 * M_PI);

	double A = mutual_firction_alpha(temp_bath) * (KAPPA*0.0001) * pow(gamma(temp_bath), 2) * pow(total_density(temp_bath) * v_n / superfluid_density(temp_bath), 3) / entropy(temp_bath);

	double T = temp_bath + A / 5 * pow(x, -5); // To change into meters.

	return T;
}