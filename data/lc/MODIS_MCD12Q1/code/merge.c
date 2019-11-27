/* *****************************************************************
 * Author: Qingsong Sun, Crystal Schaaf
 *
 * Merge MODIS products in tiles to a global data
 * Input: MODIS data in HDF format
 * Output: plain binary format
 *
 * *****************************************************************/

#include <stdio.h>
#include <assert.h>
//#include "hdf_api.c"
#include "hdf.h"
#include "mfhdf.h"

#define MAX_LEN 200

int g_global_rows;
int g_global_cols;
int g_tile_rows;
int g_tile_cols;
int32 g_data_type;
int16 g_fill_value;
char g_sds_name[MAX_LEN];
void *g_buf;
void *g_buf_tile;

int get_info(char *hdf, int dset_i);
int mosaic_tile(char *hdf, int dset_i, int h);
int process(int dset_i, int n_files, char *file_list[], char *out_name);
int get_tile_num(char *hdf, int *h, int *v);
void free_mem();
int alloc_mem();

int main(int argc, char *argv[])
{
	char *out_file = argv[1];
	int dset_i = atoi(argv[2]);
	int n_files; 
	char **file_list; 

	if(dset_i < 0){
		n_files = atoi(argv[4]);
		file_list = &argv[5];
		char *sds_name = argv[3];
		assert(n_files > 0);
		int32 id = SDstart(file_list[0], DFACC_READ);
		assert(id != FAIL);
		int32 index = SDnametoindex(id, sds_name);
		assert(index >= 0);
		dset_i = index;
	}
	else{
		n_files = atoi(argv[3]);
		file_list = &argv[4];
	}

	if(0 != process(dset_i, n_files, file_list, out_file)){
		printf("FAILD.\n");
		return -1;
	}

	return 0;
}

int select_tile(char *file_list[], int n_files, int h, int v)
{
	int index = -1;
	int _h, _v;

	int i;
	for(i=0; i<n_files; i++){
		if( 0 == get_tile_num(file_list[i], &_h, &_v)){
			if(_h == h && _v == v){
				index = i;
				break;
			}
		}
	}

	return index;
}

int process(int dset_i, int n_files, char *file_list[], char *out_name)
{
	if(0 != get_info(file_list[0], dset_i)){
		printf("Get file info failed. %s\n", file_list[0]);
		return -1;
	}

	alloc_mem();

	int32 id_out;
	int32 sds_out;

	int range[2] = {0,g_fill_value};
	int32 status;
/*
	sds_out = create_hdfeos(out_name, g_sds_name, g_sds_name, g_data_type, &g_fill_value, range,
													1.0, g_global_rows, g_global_cols, g_tile_rows, g_global_cols, &id_out);

	if(sds_out <= 0){
		printf("Create output failed. %s\n", out_name);
		free_mem();
		return -1;
	}
*/

	FILE *fp = fopen(out_name,"wb");
	if(fp == NULL){
		printf("Create output failed. %s\n", out_name);
		free_mem();
		return -1;
	}
	
	printf("Output created. %s\n", out_name);

	int h, v, i,j;

	
	int index;

	int32 start[2] = {0,0};
	int32 edge[2] = {g_tile_rows, g_global_cols};

	size_t size;
	size_t count = g_tile_rows * g_global_cols;

	for(v=0; v<18; v++){
		if(g_data_type == DFNT_UINT8){
			size = 1;
			uint8 *buf = (uint8 *)g_buf;			
			for(i=0; i<g_global_cols*g_tile_rows; i++){
				buf[i] = g_fill_value;
				//buf[i] = 0;	// ocean water
			}
		}
		else if (g_data_type == DFNT_INT16){
			size = 2;
			int16 *buf = (int16 *)g_buf;
			for(i=0; i<g_global_cols*g_tile_rows; i++){
				buf[i] = g_fill_value;
				//buf[i] = 0;	// ocean water
			}
		}
		else if (g_data_type == DFNT_UINT16){
			size = 2;
			uint16 *buf = (uint16 *)g_buf;
			for(i=0; i<g_global_cols*g_tile_rows; i++){
				buf[i] = g_fill_value;
				//buf[i] = 0;	// ocean water
			}
		}
		else{
			printf("Currently only support uint8 and int16!\n");
			return -1;
		}

		for(h=0; h<36; h++){
			index = select_tile(file_list, n_files, h, v);
			if(-1 != index){
				if(0 != mosaic_tile(file_list[index], dset_i, h)){
					printf("Write tile failed. %s\n", file_list[index]);
					free_mem();
					return -1;
				}
			}
		}
/*
		start[0] = v * g_tile_rows;
		printf("Writing %d\n", start[0]);
		status = SDwritedata(sds_out, start, NULL, edge, g_buf);
		if(status != SUCCEED){
			printf("Write Failed, v=%d\n", v);
			return -1;
		}
*/

		
		printf("Writing size=%d, count=%d v=%d\n", size, count, v);
		if(count != fwrite(g_buf, size, count, fp)){
			printf("Write Failed, v=%d\n", v);
			return -1;
		}
	}

	fclose(fp);

	printf("Succeed.\n");
	free_mem();
	return 0;
}

int get_info(char *hdf, int dset_i)
{
	char sds_name[50];
	int32 rank;
	int32 dimsizes[4];
	int32 dtype;
	int32 num_attr;

	int32 status;

	int32 id;
	int32 sds;

	id = SDstart(hdf, DFACC_READ);
	if(id <= 0){
		return -1;
	}

	sds = SDselect(id, dset_i);
	if(sds <= 0){
		SDend(id);
		return -1;
	}

	status = SDgetinfo(sds, sds_name, &rank, dimsizes, &dtype, &num_attr);
	//if(status == -1 || rank != 2){
	if(status == -1){
		SDendaccess(sds);
		SDend(id);
		return -1;
	}	

	status = SDgetfillvalue(sds, &g_fill_value);
	if(status != SUCCEED){
		g_fill_value = 255;
	}

	strcpy(g_sds_name, sds_name);
	g_global_rows = dimsizes[0] * 18;
	g_global_cols = dimsizes[1] * 36;
	g_tile_rows = dimsizes[0];
	g_tile_cols = dimsizes[1];
	g_data_type = dtype;

	SDendaccess(sds);
	SDend(id);
	return 0;
}

int alloc_mem()
{
	int type_size = 0;

	switch(g_data_type){
		case DFNT_INT8:
			;
		case DFNT_UINT8:
			type_size = 1;
			break;
		case DFNT_INT16:
			;
		case DFNT_UINT16:
			type_size = 2;
			break;
		case DFNT_INT32:
			;
		case DFNT_UINT32:
			;
		case DFNT_FLOAT32:
			type_size = 4;
			break;
		case DFNT_FLOAT64:
			type_size = 8;
			break;
		default:
			type_size = 0;
			break;
	}

	if(type_size == 0){
		printf("Error! Data type = %d\n", g_data_type);
		return -1;
	}

	g_buf = (void *)malloc(sizeof(char)*type_size*g_global_cols*g_tile_rows);
	if(g_buf == NULL){
		printf("Error! alloc failed.\n");
		return -1;
	}
	
	g_buf_tile = (void *)malloc(sizeof(char)*type_size*g_tile_cols*g_tile_rows);
	if(g_buf_tile == NULL){
		printf("Error! alloc failed.\n");
		return -1;
	}

	return 0;
}

void free_mem()
{
	if(g_buf != NULL)
		free(g_buf);
	if(g_buf_tile != NULL)
		free(g_buf_tile);
}

int is_number(char c)
{
	if(c >= 48 && c <= 57)
		return 1;
	else
		return 0;
}

int get_tile_num(char *hdf, int *h, int *v)
{
	// fine '/'
	int len = strlen(hdf);
	if(len < 6)
		return -1;

	int i;
	int pos = -1;

	for(i=len-1; i>=0; i--){
		if(hdf[i] == '/'){
			pos = i;
			break;
		}
	}

	int start;

	if(pos == -1){
		start = 0;
	}
	else{
		start = pos+1;
	}
		
	for(i=start; i<len; i++){
		if(hdf[i] == 'h'){
			if(hdf[i+3] == 'v'){
				if(is_number(hdf[i+1]) && is_number(hdf[i+2])){
					if(is_number(hdf[i+4]) && is_number(hdf[i+5])){
						*h = (hdf[i+1]-48)*10+(hdf[i+2]-48);
						*v = (hdf[i+4]-48)*10+(hdf[i+5]-48);
						return 0;
					}	
				}
			}
		}
	}

	return 1;
}

int mosaic_tile(char *hdf, int dset_i, int h)
{
	printf("Reading %s...\n", hdf);

	int start_col = h * g_tile_cols;

	if(start_col >= g_global_cols){
		printf("Write position error.Col=%d\n", start_col);
		return -1;
	}

	int32 sds_tile;
	int32 id_tile;

	id_tile = SDstart(hdf, DFACC_READ);
	if(id_tile <= 0){
		printf("Open failed. %s\n", hdf);
		return -1;
	}

	sds_tile = SDselect(id_tile, dset_i);
	if(sds_tile <= 0){
		printf("Select failed. %s, %d\n", hdf, dset_i);
		return -1;
	}

	//int32 start[2] = {0, 0};
	//int32 edge[2] = {g_tile_rows, g_tile_cols};
	// for MISR
	// R
	int32 start[4] = {0, 0, 2, 4};
	// G
	//int32 start[4] = {0, 0, 1, 4};
	// B
	//int32 start[4] = {0, 0, 0, 4};

	int32 edge[4] = {g_tile_rows, g_tile_cols, 1, 1};

	int32 status;
	
	status = SDreaddata(sds_tile, start, NULL, edge, g_buf_tile);
	if(status != SUCCEED){
		printf("Read failed. %s\n", hdf);
		SDendaccess(sds_tile);
		SDend(id_tile);
		return -1;
	}
	
	int row, col;

	for(row=0; row<g_tile_rows; row++){
		for(col=0; col<g_tile_cols; col++){
			if(g_data_type == DFNT_UINT8){
				((uint8 *)g_buf)[row*g_global_cols+col+start_col] = ((uint8 *)g_buf_tile)[row*g_tile_cols+col];
			}else if (g_data_type == DFNT_INT16){
				((int16 *)g_buf)[row*g_global_cols+col+start_col] = ((int16 *)g_buf_tile)[row*g_tile_cols+col];
			}else if (g_data_type == DFNT_UINT16){
				((uint16 *)g_buf)[row*g_global_cols+col+start_col] = ((uint16 *)g_buf_tile)[row*g_tile_cols+col];
			}
		}
	}

	SDendaccess(sds_tile);
	SDend(id_tile);

	return 0;
}
