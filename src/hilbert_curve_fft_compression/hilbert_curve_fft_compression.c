#include <hilbert_curve_fft_compression.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>



#define MIN_AMPLITUDE 5



static uint16_t _sort_partition(epicycle_t* a,uint16_t i,uint16_t j){
	float p=(a+j)->r;
	uint16_t o=i;
	while (i<j){
		if ((a+i)->r>p){
			epicycle_t tmp=*(a+o);
			*(a+o)=*(a+i);
			*(a+i)=tmp;
			o++;
		}
		i++;
	}
	epicycle_t tmp=*(a+o);
	*(a+o)=*(a+i);
	*(a+i)=tmp;
	return o;
}



static void _sort(epicycle_t* a,uint16_t i,uint16_t j){
	if (i<j){
		uint16_t k=_sort_partition(a,i,j);
		if (k){
			_sort(a,i,k-1);
		}
		_sort(a,k+1,j);
	}
}



void compress_matrix(uint8_t* bf,compressed_data_t* o){
	epicycle_t e[MATRIX_SIZE*MATRIX_SIZE];
	uint16_t i=0;
	float a=0;
	for (uint16_t j=0;j<MATRIX_SIZE*MATRIX_SIZE;j++){
		float re=0;
		float im=0;
		for (uint16_t k=0;k<MATRIX_SIZE*MATRIX_SIZE;k++){
			uint16_t x=0;
			uint16_t y=0;
			uint16_t t=k;
			for (uint8_t s=1;s<MATRIX_SIZE;s<<=1){
				uint8_t ry=t&1;
				t>>=1;
				uint8_t rx=t&1;
				t>>=1;
				if (ry^rx){
					y+=s;
				}
				else{
					if (rx){
						x=s-x-1;
						y=s-y-1;
					}
					uint16_t t=x;
					x=y;
					y=t;
				}
				x+=s*rx;
			}
			re+=(*(bf+y*MATRIX_SIZE+x))*cosf(a*k);
			im+=(*(bf+y*MATRIX_SIZE+x))*sinf(a*k);
		}
		a+=6.283185307179586f/(MATRIX_SIZE*MATRIX_SIZE);
		if (!j){
			o->off=sqrtf(re*re+im*im)/(MATRIX_SIZE*MATRIX_SIZE);
			continue;
		}
		re/=MATRIX_SIZE*MATRIX_SIZE;
		im/=MATRIX_SIZE*MATRIX_SIZE;
		e[i].r=re*re+im*im;
		if (e[i].r>=MIN_AMPLITUDE*MIN_AMPLITUDE){
			e[i].r=sqrtf(e[i].r);
			e[i].s=atan2f(-im,re);
			e[i].f=j;
			i++;
		}
	}
	if (i){
		uint16_t p=_sort_partition(e,0,i-1);
		if (p){
			_sort(e,0,p-1);
		}
		_sort(e,p+1,i-1);
	}
	o->dt=malloc(i*sizeof(epicycle_t));
	memcpy(o->dt,e,i*sizeof(epicycle_t));
	o->sz=i;
}



void decompress_matrix(compressed_data_t* c_dt,uint8_t* o){
	float a=0;
	for (uint16_t i=0;i<MATRIX_SIZE*MATRIX_SIZE;i++){
		uint16_t x=0;
		uint16_t y=0;
		uint16_t t=i;
		for (uint8_t s=1;s<MATRIX_SIZE;s<<=1){
			uint8_t ry=t&1;
			t>>=1;
			uint8_t rx=t&1;
			t>>=1;
			if (ry^rx){
				y+=s;
			}
			else{
				if (rx){
					x=s-x-1;
					y=s-y-1;
				}
				uint16_t t=x;
				x=y;
				y=t;
			}
			x+=s*rx;
		}
		float v=c_dt->off;
		for (uint16_t j=0;j<c_dt->sz;j++){
			v+=(c_dt->dt+j)->r*cosf((c_dt->dt+j)->f*a+(c_dt->dt+j)->s);
		}
		*(o+y*MATRIX_SIZE+x)=(v<0?0:(v>255?255:(uint8_t)roundf(v)));
		a+=6.283185307179586f/(MATRIX_SIZE*MATRIX_SIZE);
	}
}
