#include <stdio.h>
#include <memory.h>
#include "DFTLib.h"

#define NCOLS    256
#define NROWS    256
#define NPIXELS  (NCOLS * NROWS)

typedef unsigned char byte;

// 이미지 파일을 읽는다.
void read_image(byte* buf, const char* fname)
{
	FILE* fp = fopen(fname, "rb");

	if(fp == NULL)
	{
		fprintf(stderr, "[%s] file open error!\n", fname);
		exit(-1);
	}

	int nread = fread(buf, 1, NPIXELS, fp);
	if(nread != NPIXELS)
	{
		fprintf(stderr, "file read error! [%d/%d]\n", nread, NPIXELS);
		exit(-1);
	}

	printf("*image read success! : %s \n", fname);
}

// 이미지 파일을 만든다.
void write_image(const byte* buf, const char* fname)
{
	FILE* fp = fopen(fname, "wb");
	if(fp == NULL)
	{
		fprintf(stderr, "[%s] file open error!\n", fname);
		exit(-1);
	}

	int nwrite = fwrite(buf, 1, NPIXELS, fp);
	if(nwrite != NPIXELS)
	{
		fprintf(stderr, "file write error! [%d/%d]\n", nwrite, NPIXELS);
		exit(-1);
	}

	printf("*image write success! : %s \n", fname);
}

// Lena 이미지를 읽고 double 형의 데이터로 반환
void load_lena_image(double* dr, double* di)
{
	byte image_input[NPIXELS] = {0,};
	read_image(image_input, "lena256.raw");

	// 원본 Pixel data는 실수 영역만 존재하고 허수는 모두 0
	for(int i=0; i<NPIXELS; i++)
	{
		dr[i] = image_input[i];
		di[i] = 0.0L;
	}
}

void write_double_image(double* d, const char* fname)
{
	byte image_output[NPIXELS] = {0,};
	for(int i=0; i<NPIXELS; i++)
	{
		image_output[i] = (byte)d[i];
	}

	write_image(image_output, fname);
}

void dft2d(double* dr, double* di)
{
	// 2 차원 DFT
	for(int j=0; j<NROWS; j++)
	{
		dft(NCOLS, dr + j*NCOLS, di + j*NCOLS);
	}

	for(int i=0; i<NCOLS; i++)
	{
		dft_y(NROWS, NCOLS, dr + i, di + i);
	}

	printf("*2D-DFT ... ok!\n");
}

void inverse_dft2d(double* dr, double* di)
{
	for(int i=0; i<NCOLS; i++)
	{
		inverse_dft_y(NROWS, NCOLS, dr + i, di + i);
	}
	
	for(int j=0; j<NROWS; j++)
	{
		inverse_dft(NCOLS, dr + j*NCOLS, di + j*NCOLS);
	}

	printf("*Inverse 2D-DFT ... ok!\n");
}

void low_pass_filtering(double* dr, double* di, int radius)
{
	double cj = NROWS/2; // 중심 y축
	double ci = NCOLS/2; // 중심 x축

	for(int j=0; j<NROWS; j++)
	{
		for(int i=0; i<NCOLS; i++)
		{
			if(sqrt((long double)((j-cj)*(j-cj) + (i-ci)*(i-ci))) > (long double)radius)
			{
				dr[j*NCOLS + i] = 0;
				di[j*NCOLS + i] = 0;
			}
		}
	}

	printf("*%d pixel low-pass filtering ... ok!\n", radius);
}

void main()
{
	double* dr = (double *)malloc(NPIXELS * sizeof(double));
	double* di = (double *)malloc(NPIXELS * sizeof(double));

	printf("--------------------------------------------------------------------------\n");
	printf("- 원본 -> 2D-DFT -> Inverse 2D-DFT\n");
	printf("--------------------------------------------------------------------------\n");

	load_lena_image(dr, di);

	dft2d(dr, di);

	write_double_image(dr, "lena256_dft2d_real.raw");
	write_double_image(di, "lena256_dft2d_imag.raw");

	inverse_dft2d(dr, di);

	write_double_image(dr, "lena256_inverse_dft2d.raw");

	printf("--------------------------------------------------------------------------\n");
	printf("- 원본 -> 2D-DFT -> 10 pixel low-pass filtering -> Inverse 2D-DFT\n");
	printf("--------------------------------------------------------------------------\n");

	load_lena_image(dr, di);

	dft2d(dr, di);

	low_pass_filtering(dr, di, 10);

	inverse_dft2d(dr, di);

	write_double_image(dr, "lena256_dft2d_10px_filtered.raw");

	printf("--------------------------------------------------------------------------\n");
	printf("- 원본 -> 2D-DFT -> 25 pixel low-pass filtering -> Inverse 2D-DFT\n");
	printf("--------------------------------------------------------------------------\n");

	load_lena_image(dr, di);

	dft2d(dr, di);

	low_pass_filtering(dr, di, 25);

	inverse_dft2d(dr, di);

	write_double_image(dr, "lena256_dft2d_25px_filtered.raw");

	printf("--------------------------------------------------------------------------\n");
	printf("- 원본 -> 2D-DFT -> 50 pixel low-pass filtering -> Inverse 2D-DFT\n");
	printf("--------------------------------------------------------------------------\n");

	load_lena_image(dr, di);

	dft2d(dr, di);

	low_pass_filtering(dr, di, 50);

	inverse_dft2d(dr, di);

	write_double_image(dr, "lena256_dft2d_50px_filtered.raw");

	printf("--------------------------------------------------------------------------\n");
	printf("- 원본 -> 2D-DFT -> 100 pixel low-pass filtering -> Inverse 2D-DFT\n");
	printf("--------------------------------------------------------------------------\n");

	load_lena_image(dr, di);

	dft2d(dr, di);

	low_pass_filtering(dr, di, 100);

	inverse_dft2d(dr, di);

	write_double_image(dr, "lena256_dft2d_100px_filtered.raw");

	free(dr);
	free(di);
}