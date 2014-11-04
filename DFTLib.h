#include <math.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

// Take, eat, this is my code which is hacked for you...
// (In case you didn't grasp that one consider this public domain)
// nolandda Thu Mar 11 03:49:18 EST 2004

/*   Discrete Fourier Transform */
int dft(long int length, double real_sample[], double imag_sample[]);
int dft_y(long int length, long int col_len, double real_sample[], double imag_sample[]);

/*   Inverse Discrete Fourier Transform */
int inverse_dft(long int length, double real_sample[], double imag_sample[]);
int inverse_dft_y(long int length, long int col_len, double real_sample[], double imag_sample[]);


int test();