#ifndef __HILBERT_CURVE_FFT_COMPRESSION_H__
#define __HILBERT_CURVE_FFT_COMPRESSION_H__ 1
#include <stdint.h>



#define MATRIX_SIZE 4



typedef struct __EPICYCLE{
	float r;
	float s;
	uint16_t f;
} epicycle_t;



typedef struct __COMPRESSED_DATA{
	epicycle_t* dt;
	uint32_t sz;
	float off;
} compressed_data_t;



void compress_matrix(uint8_t* bf,compressed_data_t* o);



void decompress_matrix(compressed_data_t* c_dt,uint8_t* o);



#endif
