#ifdef _MSC_VER
#include <windows.h>
#endif
#include <hilbert_curve_fft_compression.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>



void print_matrix(uint8_t* m){
	uint16_t i=0;
	while (i<MATRIX_SIZE*MATRIX_SIZE){
		printf("\x1b[48;2;%u;%u;%um  ",m[i],m[i],m[i]);
		i++;
		if (!(i&(MATRIX_SIZE-1))){
			printf("\x1b[0m\n");
		}
	}
}



int main(int argc,const char** argv){
#ifdef _MSC_VER
	SetConsoleOutputCP(CP_UTF8);
	SetConsoleMode(GetStdHandle(-11),7);
#endif
	uint8_t m[MATRIX_SIZE*MATRIX_SIZE]={
		255,255,255,255,
		128,64,23,125,
		102,75,163,244,
		255,67,218,255
	};
	compressed_data_t dt;
	print_matrix(m);
	compress_matrix(m,&dt);
	printf("Offset: %f\n",dt.off);
	for (uint16_t i=0;i<dt.sz;i++){
		printf("[%u]: %u, %f, %f\n",i,(dt.dt+i)->f,(dt.dt+i)->r,(dt.dt+i)->s);
	}
	decompress_matrix(&dt,m);
	free(dt.dt);
	print_matrix(m);
	return 0;
}
