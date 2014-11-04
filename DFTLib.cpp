/*                               -*- Mode: C -*- 
* Filename: dft.c
* Copyright (C) Dan Noland 2003
* Author: Dan Noland
* Created: Thu Mar 11 03:19:19 2004
*           By: Dan Noland
* Last-Updated: Thu Mar 11 03:51:54 2004
*     Update #: 12
* Status: 
*/

//#include "stdafx.h"
#include <stdio.h>
#include "DFTLib.h"

#define NEW

// Take, eat, this is my code which is hacked for you...
// (In case you didn't grasp that one consider this public domain)
// nolandda Thu Mar 11 03:49:18 EST 2004

/*
Discrete Fourier Transform
*/

int dft(long int length, double real_sample[], double imag_sample[])
{
	long int i, j;
	double arg;
	double cosarg,sinarg;
	double *temp_real=NULL,*temp_imag=NULL;

	temp_real = (double *)malloc(length *  sizeof(double));
	temp_imag = (double *)malloc(length *  sizeof(double));
	if (temp_real == NULL || temp_imag == NULL)
	{
		return(FALSE);
	}

	for(i=0; i<length; i+=1) 
	{
		temp_real[i] = 0;
		temp_imag[i] = 0;
		arg = -1.0 * 2.0 * 3.141592654 * (double)i / (double)length;
		for(j=0; j<length; j+=1) 
		{
			cosarg = cos(j * arg);
			sinarg = sin(j * arg);
			temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
			temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
		}
	}

	/* Copy the data back */
	for (i=0; i<length; i+=1) 
	{
#ifndef NEW
		real_sample[i] = temp_real[i];
		imag_sample[i] = temp_imag[i];
#else
		real_sample[i] = temp_real[(i+length/2)%(length)];
		imag_sample[i] = temp_imag[(i+length/2)%(length)];
#endif
	}

	free(temp_real);
	free(temp_imag);
	return(TRUE);
}

// ���ι��� dft
int dft_y(long int length, long int col_len, double real_sample[], double imag_sample[])
{
	long int i, j;
	double arg;
	double cosarg,sinarg;
	double *temp_real=NULL,*temp_imag=NULL;

	temp_real = (double *)malloc(length *  sizeof(double));
	temp_imag = (double *)malloc(length *  sizeof(double));
	if (temp_real == NULL || temp_imag == NULL)
	{
		return(FALSE);
	}

	for(i=0; i<length; i+=1) 
	{
		temp_real[i] = 0;
		temp_imag[i] = 0;
		arg = -1.0 * 2.0 * 3.141592654 * (double)i / (double)length;
		for(j=0; j<length; j+=1) 
		{
			cosarg = cos(j * arg);
			sinarg = sin(j * arg);
			temp_real[i] += (real_sample[j*col_len] * cosarg - imag_sample[j*col_len] * sinarg);
			temp_imag[i] += (real_sample[j*col_len] * sinarg + imag_sample[j*col_len] * cosarg);
		}
	}

	/* Copy the data back */
	for (i=0; i<length; i+=1) 
	{
#ifndef NEW
		real_sample[i*col_len] = temp_real[i];
		imag_sample[i*col_len] = temp_imag[i];
#else
		real_sample[i*col_len] = temp_real[(i+length/2)%(length)];
		imag_sample[i*col_len] = temp_imag[(i+length/2)%(length)];
#endif
	}

	free(temp_real);
	free(temp_imag);
	return(TRUE);
}


/*
Inverse Discrete Fourier Transform
*/

int inverse_dft(long int length, double real_sample[], double imag_sample[])
{
	long int i, j;
	double arg;
	double cosarg,sinarg;
	double *temp_real=NULL,*temp_imag=NULL;

	temp_real = (double *)malloc(length * sizeof(double));
	temp_imag = (double *)malloc(length * sizeof(double));
	if (temp_real == NULL || temp_imag == NULL)
	{
		return(FALSE);
	}

#ifdef NEW
	for (i=0; i<length; i+=1) 
	{
		temp_real[i] = real_sample[i];
		temp_imag[i] = imag_sample[i];
	}
	for (i=0; i<length; i+=1) 
	{
		real_sample[i] = temp_real[(i+length/2)%length];
		imag_sample[i] = temp_imag[(i+length/2)%length];
	}
#endif

	for(i=0; i<length; i+=1) 
	{
		temp_real[i] = 0;
		temp_imag[i] = 0;
		arg = 2.0 * 3.141592654 * (double)i / (double)length;
		for(j=0; j<length; j+=1) 
		{
			cosarg = cos(j * arg);
			sinarg = sin(j * arg);
			temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
			temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
		}
	}

	/* Copy the data back */
	for (i=0; i<length; i+=1) 
	{
		real_sample[i] = temp_real[i] / (double)length;
		imag_sample[i] = temp_imag[i] / (double)length;
	}

	free(temp_real);
	free(temp_imag);
	return(TRUE);
}

// ���ι��� inverse dft
int inverse_dft_y(long int length, long int col_len, double real_sample[], double imag_sample[])
{
	long int i, j;
	double arg;
	double cosarg,sinarg;
	double *temp_real=NULL,*temp_imag=NULL;

	temp_real = (double *)malloc(length * sizeof(double));
	temp_imag = (double *)malloc(length * sizeof(double));
	if (temp_real == NULL || temp_imag == NULL)
	{
		return(FALSE);
	}

#ifdef NEW
	for (i=0; i<length; i+=1) 
	{
		temp_real[i] = real_sample[i*col_len];
		temp_imag[i] = imag_sample[i*col_len];
	}
	for (i=0; i<length; i+=1) 
	{
		real_sample[i*col_len] = temp_real[(i+length/2)%length];
		imag_sample[i*col_len] = temp_imag[(i+length/2)%length];
	}
#endif

	for(i=0; i<length; i+=1) 
	{
		temp_real[i] = 0;
		temp_imag[i] = 0;
		arg = 2.0 * 3.141592654 * (double)i / (double)length;
		for(j=0; j<length; j+=1) 
		{
			cosarg = cos(j * arg);
			sinarg = sin(j * arg);
			temp_real[i] += (real_sample[j*col_len] * cosarg - imag_sample[j*col_len] * sinarg);
			temp_imag[i] += (real_sample[j*col_len] * sinarg + imag_sample[j*col_len] * cosarg);
		}
	}

	/* Copy the data back */
	for (i=0; i<length; i+=1) 
	{
		real_sample[i*col_len] = temp_real[i] / (double)length;
		imag_sample[i*col_len] = temp_imag[i] / (double)length;
	}

	free(temp_real);
	free(temp_imag);
	return(TRUE);
}

int test()
{
	int i;
	double* dr = (double *)malloc(100 * sizeof(double));
	double* di = (double *)malloc(100 * sizeof(double));

	/* ���� Pixel data�� �Ǽ� ������ �����ϰ� ����� ��� 0 */
	for(i=0; i<100; i++)
	{
		dr[i] = i;
		di[i] = 0;
	}

	/* DFT ���� ����� Ȯ�� ���� */
	dft(100, dr, di);
	for(i=0; i<100; i++)
	{
		printf("[%3d] %10.3f %10.3f\n", i, dr[i], di[i]);
	}

	/* Inverse DFT �Ŀ� ��ȣ�� ���������� �����Ǵ� ���� Ȯ�� */
	inverse_dft(100, dr, di);
	printf("\n\n reconstructed....\n");
	for(i=0; i<100; i++)
	{
		printf("[%3d] %10.3f %10.3f\n", i, dr[i], di[i]);
	}

	return 1;
}