/**
	Header for frequency filter smoothing high spatial frequencies.
	@file vortex_smooth.h
*/

#include "vortex_smooth.h"

#include <stdio.h>
#include <assert.h>
#include "math.h"
#include <fftw3.h>

#include "vec3_math.h"
#include "tangle.h"
#include "vortex_constants.h"
#include "boundaries.h"

/**
	Find out if a vertox containing given point is open or closed. If the vortex is open returns starting point, if closed returns -1.
	@param tangle: Tangle structure holding all information about vortices.
	@param vortex_point: Point of te vortex we want to find information about.
	@return Returns starting point if the vortex is open, returns -1 if closed.
*/
int is_vortex_open(struct tangle_state *tangle, int vortex_point)
{
    int i = vortex_point;
	// Loop throught the vortex, while there is where to go or we haven't got to the starting point.
    while(tangle->connections[i].reverse != vortex_point && tangle->connections[i].reverse != -1)
    {
        i = tangle->connections[i].reverse;
    }

	// Returns result.
    if (tangle->connections[i].reverse == vortex_point) return -1;
    if (tangle->connections[i].reverse == -1) return i;

	printf("\nERROR in is_vortex_open.\n");
	return EXIT_FAILURE;
}

/**
	Find out length of the open vortex and its ending point.
	@param tangle: Tangle structure holding all information about vortices.
	@param starting_point: Starting point of the vortex.
	@param ending_point: Pointer to ending point of the vortex.
	@param length: Pointer to the length of the vortex.
*/
void find_ending_point(struct tangle_state *tangle, int starting_point, int *ending_point, int *length)
{
    int i = starting_point;
    *length = 1;
	// Loop through the vortex with right direction.
	if (tangle->connections[i].forward != -1)
	{
		while (tangle->connections[i].forward != -1)
		{
			i = tangle->connections[i].forward;
			(*length)++;
		}
	}
	else
	{
		while(tangle->connections[i].reverse != -1)
		{
			i = tangle->connections[i].reverse;
			(*length)++;
		}
	}
    
    *ending_point = i;
}

/**
	Find out length of the closed vortex.
	@param tangle: Tangle structure holding all information about vortices.
	@param starting_point: Starting point of the vortex.
	@returns: Length of the vortex.
*/
int closed_vortex_length(struct tangle_state *tangle, int starting_point)
{
    int i = starting_point;
    int length = 1;
	// Loop through the vortex.
    while(tangle->connections[i].forward != starting_point)
    {
        i = tangle->connections[i].forward;
        length++;
    }

    return length;
}

/*
	Prints error call if the distance of node from if original position is greater then 5 * dl_max
*/
void check_distance(struct tangle_state* tangle, struct vec3 r1, struct vec3 r2) {
	struct segment seg = seg_pwrap(&r1, &r2, &tangle->box);

	if (segment_len(&seg) > 5 * global_dl_max) printf("\nERROR: Smooth moves vortex away from original position by distance greater then 5 * dl_max.\n");
}

/**
	!!! UNUSED !!!
*/
void fill_coordinate_array(double *array, int length, int side, int parity)
{
	/// BEWARE! actall length of array is length-1 + length + length-1
	/// if you want to use this function, rewrite it !!!
    // side: 0 = beggining, 1 = ending
    // parity: 0 = even, 1 = odd
    switch (10*side+parity) {
        case 0:
			for (int i = 0; i < length - 1; i++) *(array + i) = *(array + 2 * length - 3 - i);
            break;
        case 1: ;
			for (int i = 0; i < length - 1; i++) *(array + i) = - *(array + 2 * length - 3 - i);
            break;
        case 10: ;
			for (int i = 0; i < length - 1; i++) *(array + 2 * length - 1 + i) = *(array + 2 * length - 3 - i);
            break;
        case 11: ;
			for (int i = 0; i < length - 1; i++) *(array + 2 * length - 1 + i) = - *(array + 2 * length - 3 - i);
            break;
        default: ;
            printf("ERROR: Unknown pasrity/side in fill_coordinate_array. (WORTEX_SMOOTH.C)");
    }
}

/**
	Smooth the periodic discrete function (array of positions) using frequency filter smoothing high spatial frequencies.
	@param in: Periodic discrete function (array) to smooth.
	@param N: Length of the array *in.
*/
void low_pass_filter(double *in, int N)
{
	// Defines array to fill with result.
	double *out;
	out = (double*)fftw_malloc(sizeof(double) * N);

	// Define variables needed for fftw. Execute forward, real to real, 1 dimensional, fast fourier transform.
	fftw_plan forward_plan, backward_plan;
	forward_plan = fftw_plan_r2r_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(forward_plan);

	// Filter high frequencies.
	int smooth_number = freq_to_cutoff * N / 2; // This is number of cutted frequencies. Don't forget that there is N/2 frequencies in total, not N! Hence 3*N/8 means 3/4 frequencies.
	if (N > 2 * smooth_number+1) {
		// i=0 constant, i=1 ... N/2 cosine (positive means +cos), i=N/2 ... N-1 sines (positive means -sin)
		if (N % 2 == 0) {
			// Odd numb. of frequencies + 1 constant.
			for (int i = N/2-smooth_number; i < N/2+1+smooth_number; i++) {
                *(out+i) = 0;
			}
		} else {
			// Even numb. of frequencies + 1 constant.
			for (int i = N/2-smooth_number; i < N/2+smooth_number; i++) {
				*(out+i) = 0;
			}
		}
	}

	// Execute backward, real to real, 1 dimensional, fast fourier transform.
	backward_plan = fftw_plan_r2r_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(backward_plan);

	// Normalize result.
	for (int i = 0; i < N; i++) {
		*(in+i) = *(in+i)/N;
	}

	// Free fftw stuff.
	fftw_destroy_plan(forward_plan);
	fftw_destroy_plan(backward_plan);
	fftw_free(out);
}

///// TODO: CHECK THE DISTANCE FROM ORIGINAL POSITION, IS IT MORE THEN 2*dl_max?

/**
	Frequency filter smoothing high spatial frequencies.
	@param tangle: Tangle structure holding all information about vortices we want to smooth.
*/
void high_freq_cutoff(struct tangle_state *tangle)
{
    // Prepare array to remember visited points.
    int visited[tangle->N];
    for (int i = 0; i < tangle->N; i++) {
        visited[i] = 0;
    }

    // Loop through vortices.
    for (int i = 0; i < tangle->N; i++) {
        // Vortex is visited or is empty => continue to another point.
        if (visited[i] == 1) continue;
        if (tangle->status[i].status == EMPTY) {
            visited[i] = 1;
            continue;
        }

        // Not yet visited => this loop is not smooth yet! Find if it is an open or a closed vortex, find starting point.
        // is_vortex_open returns -1 for closed vortex, returns starting point for open vortex
        int starting_point = is_vortex_open(tangle, i);

		// Case of the closed vortex.
        if (starting_point == -1)
		{
            starting_point = i;

			// Finds length of the closed vortex.
            int length = closed_vortex_length(tangle, starting_point);

            // Prepare coordinate arrays.
            double *x = (double*)malloc((length)*sizeof(double));
            double *y = (double*)malloc((length)*sizeof(double));
            double *z = (double*)malloc((length)*sizeof(double));

            // Fill the 0-th element.
            int k = starting_point;
            x[0] = tangle->vnodes[k].p[0];
            y[0] = tangle->vnodes[k].p[1];
            z[0] = tangle->vnodes[k].p[2];
			// Mark the 0-th element as visited.
            visited[k] = 1;
			
			// Move forward.
            k = tangle->connections[k].forward;

            // Fill rest of the coordinate arrays using pwrap (it is periodic so we can fill whole array).
            struct segment seg;
            struct vec3 k_old;
			// Actual loop.
            for (int j = 1; j < length; j++) {
                k_old = vec3(x[j-1], y[j-1], z[j-1]);
                seg = seg_pwrap(&k_old, tangle->vnodes + k, &tangle->box);
                x[j] = seg.r2.p[0];
                y[j] = seg.r2.p[1];
                z[j] = seg.r2.p[2];
				// Mark elements as visited.
                visited[k] = 1;
				// Move forward.
                k = tangle->connections[k].forward;
            }

			// Look if the vortex is infinitely periodic.
			struct vec3 start_to_end = vec3(x[0]-x[length-1], y[0]-y[length-1], z[0]-z[length-1]);
			if (vec3_len(&start_to_end) > (tangle->box.top_right_front.p[0] - tangle->box.bottom_left_back.p[0]) / 2)
			{
				// Increase the size of the array.
				x = (double*)realloc(x, (2*length-2) * sizeof(double));
				y = (double*)realloc(y, (2*length-2) * sizeof(double));
				z = (double*)realloc(z, (2*length-2) * sizeof(double));
				// Subtract coordinates from direct line between points.
				// Linear coefficients and constant.
				double x_coef = (x[length - 1] - x[0]) / (length - 1);
				double y_coef = (y[length - 1] - y[0]) / (length - 1);
				double z_coef = (z[length - 1] - z[0]) / (length - 1);
				double x_0 = x[0];
				double y_0 = y[0];
				double z_0 = z[0];
				//Bboundary points with respect to the line are zero.
				x[length - 1] = 0;
				y[length - 1] = 0;
				z[length - 1] = 0;
				x[0] = 0;
				y[0] = 0;
				z[0] = 0;
				// Actual loop for subtracting.
				for (int j = 1; j < length - 1; j++) {
					x[j] = x[j] - x_coef * j - x_0;
					y[j] = y[j] - y_coef * j - y_0;
					z[j] = z[j] - z_coef * j - z_0;
				}

				// Fill up rest of coordinate arrays to make it periodic (needed for fftw).
				for (int j = 0; j < length - 2; j++) {
					x[length + j] = -x[length - 2 - j];
					y[length + j] = -y[length - 2 - j];
					z[length + j] = -z[length - 2 - j];
				}

				// Do low pass filter.
				low_pass_filter(x, length * 2 - 2);
				low_pass_filter(y, length * 2 - 2);
				low_pass_filter(z, length * 2 - 2);

				// Redoing back subtraction from line.
				for (int j = 0; j < length; j++) {
					x[j] = x[j] + x_coef * j + x_0;
					y[j] = y[j] + y_coef * j + y_0;
					z[j] = z[j] + z_coef * j + z_0;
				}
			}
			// The vortex is classic loop.
			else {
				// Do low pass filter.
				low_pass_filter(x, length);
				low_pass_filter(y, length);
				low_pass_filter(z, length);
			}

            // Rewrite coordinates back to tangle.
            k = starting_point;
            for (int j = 0; j < length; j++) {
				check_distance(tangle, vec3(x[j], y[j], z[j]), tangle->vnodes[k]);
                tangle->vnodes[k].p[0] = x[j];
                tangle->vnodes[k].p[1] = y[j];
                tangle->vnodes[k].p[2] = z[j];
                k = tangle->connections[k].forward;
            }

            // Free arrays.
            free(x);
            free(y);
            free(z);
        }
		// Case of the open vortex.
		else
		{
			// Finds the ending point and the length of the open vortex.
			int ending_point;
            int length;
            find_ending_point(tangle, starting_point, &ending_point, &length);

			// Prepare the coordinate arrays, they are longer becouse of mirror copy of array exept it endpoints (thats why 2*length - 2, that -2 is for endpoints.)
			double* x = (double*)malloc((length * 2 - 2) * sizeof(double));
			double* y = (double*)malloc((length * 2 - 2) * sizeof(double));
			double* z = (double*)malloc((length * 2 - 2) * sizeof(double));

			// Fill the 0-th element.
			int k = starting_point;
			x[0] = tangle->vnodes[k].p[0];
			y[0] = tangle->vnodes[k].p[1];
			z[0] = tangle->vnodes[k].p[2];
			// Mark the 0-th element as visited.
			visited[k] = 1;
			// Move forward.
			k = tangle->connections[k].forward;

			// Fill rest of the points using pwrap.
			struct segment seg;
			struct vec3 k_old;
			for (int j = 1; j < length; j++) {
				k_old = vec3(x[j - 1], y[j - 1], z[j - 1]);
				seg = seg_pwrap(&k_old, tangle->vnodes + k, &tangle->box);
				x[j] = seg.r2.p[0];
				y[j] = seg.r2.p[1];
				z[j] = seg.r2.p[2];
				// Mark the k-th element as visited.
				visited[k] = 1;
				// Move forward.
				k = tangle->connections[k].forward;
			}

            // Subtract coordinates from direct line between points.
            // Linear coefficients and constant.
            double x_coef = (x[length - 1] - x[0]) / (length - 1);
            double y_coef = (y[length - 1] - y[0]) / (length - 1);
            double z_coef = (z[length - 1] - z[0]) / (length - 1);
            double x_0 = x[0];
            double y_0 = y[0];
            double z_0 = z[0];
            //Bboundary points with respect to the line are zero.
            x[length - 1] = 0;
            y[length - 1] = 0;
            z[length - 1] = 0;
            x[0] = 0;
            y[0] = 0;
            z[0] = 0;
            // Actual loop for subtracting.
            for (int j = 1; j < length-1; j++) {
                x[j] = x[j] - x_coef*j - x_0;
                y[j] = y[j] - y_coef*j - y_0;
                z[j] = z[j] - z_coef*j - z_0;
            }

            // Fill up rest of coordinate arrays to make it periodic (needed for fftw).
			for (int j = 0; j < length - 2; j++) {
				x[length + j] = -x[length - 2 - j];
				y[length + j] = -y[length - 2 - j];
				z[length + j] = -z[length - 2 - j];
			}


			/// UNUSED ///
			/*
            int start_pin_mode = tangle->status[starting_point].status;
            int end_pin_mode = tangle->status[starting_point].status;
            int start_pin_wall = tangle->status[starting_point].pin_wall;
            int end_pin_wall = tangle->status[starting_point].pin_wall;

            if (start_pin_mode == PINNED || start_pin_wall == LEFT  || start_pin_wall == RIGHT) {
                fill_coordinate_array(x, length, 0, 1);
            } else {
                fill_coordinate_array(x, length, 0, 0);
            }
            if (start_pin_mode == PINNED || start_pin_wall == FRONT  || start_pin_wall == BACK) {
                fill_coordinate_array(y, length, 0, 1);
            } else {
                fill_coordinate_array(y, length, 0, 0);
            }
            if (start_pin_mode == PINNED || start_pin_wall == UP  || start_pin_wall == DOWN) {
                fill_coordinate_array(z, length, 0, 1);
            } else {
                fill_coordinate_array(z, length, 0, 0);
            }
            if (end_pin_mode == PINNED || end_pin_wall == LEFT  || end_pin_wall == RIGHT) {
                fill_coordinate_array(x, length, 1, 1);
            } else {
                fill_coordinate_array(x, length, 1, 0);
            }
            if (end_pin_mode == PINNED || end_pin_wall == FRONT  || end_pin_wall == BACK) {
                fill_coordinate_array(y, length, 1, 1);
            } else {
                fill_coordinate_array(y, length, 1, 0);
            }
            if (end_pin_mode == PINNED || end_pin_wall == UP  || end_pin_wall == DOWN) {
                fill_coordinate_array(z, length, 1, 1);
            } else {
                fill_coordinate_array(z, length, 1, 0);
            }
			*/

			// Do low pass filter.
			low_pass_filter(x, length*2-2);
			low_pass_filter(y, length*2-2);
			low_pass_filter(z, length*2-2);

			// Redoing back subtraction from line.
			for (int j = 0; j < length; j++) {
				x[j] = x[j] + x_coef * j + x_0;
				y[j] = y[j] + y_coef * j + y_0;
				z[j] = z[j] + z_coef * j + z_0;
			}

			// Rewrite coordinates to the tangle, exept boundary points.
			k = starting_point;
			k = tangle->connections[k].forward;
			for (int j = 1; j < length - 1; j++) {
				check_distance(tangle, vec3(x[j], y[j], z[j]), tangle->vnodes[k]);
				tangle->vnodes[k].p[0] = x[j];
				tangle->vnodes[k].p[1] = y[j];
				tangle->vnodes[k].p[2] = z[j];
				k = tangle->connections[k].forward;
			}

			// Free arrays.
			free(x);
			free(y);
			free(z);

		}
    }
}

/**
	Main call function for frequency filter smoothing high spatial frequencies.
	@param tangle: Tangle structure holding all information about vortices we want to smooth.
*/
void smooth(struct tangle_state *tangle)
{
    if (use_freq_filter) high_freq_cutoff(tangle);
}
