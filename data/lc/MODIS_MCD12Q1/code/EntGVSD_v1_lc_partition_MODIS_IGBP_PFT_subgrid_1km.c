/***********************************************************************************
 * Author: Qingsong Sun, Crystal Schaaf, 2016-05-01                                *
 *                                                                                 *
 * Ent Global Vegetation Structure Data set v1 processing step.                    *
 *                                                                                 *
 * Takes the outputs from EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_500m.exe          *
 * consisting of 500 m gridded maps of 29 cover classes and organizes them into    *
 * subgrid cover fractions at 1 km.  Outputs are in NetCDF format.                 *
 * Results are used as input files for further classification into Ent PFTs and    *
 * land cover types.                                                               *
 *                                                                                 *
 * ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "netcdf.h"
#include <assert.h>
#include <string.h>

#define HDFEOS_GROUP_NAME "IGBP_PFT_PARTITION"

#define N_SDS (29)

#define COAST_LINE (32767)

#define ERRCODE 2
#define ERR(e) {if(e){printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}}

int nrows = 21600;
int ncols = 43200;

int create_nc(char *file_name, char *var_name, int *_ncid, int *_varid);

int main(int argc, char *argv[])
{
	if(argc != 3){
		printf("usage: %s partition.bin out.hdf coast.bin\n", argv[0]);
		return 0;
	}

	char *file_part = argv[1];
	char *nc_out = argv[2];

	int i, j;
	char sds_names[N_SDS][100];
	char nc_names[N_SDS][100];

	int ncids[N_SDS];
	int varids[N_SDS];
	int ret;

	for(i=0; i<N_SDS; i++){
		sprintf(sds_names[i], "PARTITION_%d", i);
		sprintf(nc_names[i], "%s.PARTITION_%02d.nc", nc_out, i);
		assert(0 == create_nc(nc_names[i], sds_names[i], &ncids[i], &varids[i]));
	}

	FILE *fp = fopen(file_part, "rb");
	if(fp == NULL){
		printf("Open %s failed.\n", file_part);
		return 0;
	}

	unsigned char buf[2][86400];
	//int32 start[2] = {0,0};
	//int32 edge[2] = {1, ncols};
	short buf_out[N_SDS][ncols];

	int row, col;
	int n_read = 86400 * 2;
	unsigned char v;

	size_t start[2] = {0, 0};
	size_t count[2] = {1, ncols};

	for(row=0; row<nrows; row++){
		start[0] = row;

		printf("\r%3d%%", 100*row/nrows);
		fflush(stdout);
		
		if(n_read != fread(buf, 1, n_read, fp)){
			printf("Read failed.\n");
			goto END;
		}
		
		for(col=0; col<ncols; col++){
			for(i=0; i<N_SDS; i++){
				buf_out[i][col] = 0;	
			}

			for(i=0; i<2; i++){
				for(j=0; j<2; j++){
					v = buf[i][j+col*2];
					if(v < N_SDS){
						buf_out[v][col] ++;	
					}
				}
			}

			for(i=0; i<N_SDS; i++){
				if(buf_out[i][col] > 0){
					buf_out[i][col] = 100 * buf_out[i][col] / 4;
				}
			}

			// draw coast line
			/*if(coast[row*1440+col] == 255){
				for(i=0; i<N_SDS; i++){
					buf_out[i][col] = COAST_LINE;
				}
			}*/


		}	// col

		for(i=0; i<N_SDS; i++){
			ret = nc_put_vara_short(ncids[i], varids[i], start, count, buf_out[i]);
			ERR(ret);
		}

	}	// row

	for(i=0; i<N_SDS; i++){
		ret = nc_close(ncids[i]);
		assert(ret == NC_NOERR);
	}

END:
	if(fp != NULL){
		fclose(fp);
		fp = NULL;
	}

	/*for(i=0; i<N_SDS; i++){
		SDendaccess(sds_ids[i]);
	}
	SDend(id);*/
	return 0;
}

int create_nc(char *file_name, char *var_name, int *_ncid, int *_varid)
{
	printf("creating %s...\n", file_name);
	int nrows = 21600;
	int ncols = 43200;
	
	int chunk_rows = 10;

	/*if(SUCCEED != create_hdfeos(hdf_out, sds_names, N_SDS, data_type, &fill_value, nrows, ncols, chunk_rows, &id, sds_ids)){
		printf("Create %s failed.\n", hdf_out);
		return 0;
	}*/

	int ncid;
	int varid;
	int lat, lon;
	int lat_dimid, lon_dimid;
	int lat_varid, lon_varid;
	int ret = nc_create(file_name, NC_64BIT_OFFSET, &ncid);
	assert(ret == NC_NOERR);
	ret = nc_def_dim(ncid, "lat", nrows, &lat_dimid);
	assert(ret == NC_NOERR);
	ret = nc_def_dim(ncid, "lon", ncols, &lon_dimid);
	assert(ret == NC_NOERR);
	ret = nc_def_var(ncid, "lat", NC_FLOAT, 1, &lat_dimid, &lat_varid);
	assert(ret == NC_NOERR);
	ret = nc_def_var(ncid, "lon", NC_FLOAT, 1, &lon_dimid, &lon_varid);
	assert(ret == NC_NOERR);

	char *unit_lat = "degrees_north";
	char *unit_lon = "degrees_east";

	ret = nc_put_att_text(ncid, lat_varid, "units", strlen(unit_lat), unit_lat);
	assert(ret == NC_NOERR);
	ret = nc_put_att_text(ncid, lon_varid, "units", strlen(unit_lon), unit_lon);
	assert(ret == NC_NOERR);

	int dimids[2] = {lat_dimid, lon_dimid};
	size_t chunks[2] = {ncols, nrows/100};
	nc_type xtype = NC_SHORT;

	short fill = 32767;
	ret = nc_def_var(ncid, var_name, xtype, 2, dimids, &varid);
	assert(ret == NC_NOERR);
	/*int shuffle = NC_SHUFFLE;
	int deflate = 1;        // This switches compression on (1) or off (0).
	int deflate_level = 8;  // This is the compression level in range 1 (less) - 9 (more).
	ret = nc_def_var_chunking(ncid, varid, 0, chunks);
	ERR(ret);
	ret = nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level);
	ERR(ret);*/

	ret = nc_put_att_short(ncid, varid, "_FillValue", NC_SHORT, 1, &fill);
	ERR(ret);
	ret = nc_enddef(ncid);
	ERR(ret);

	float lats[nrows];
	float lons[ncols];

	// global
	int i;
	for(i=0; i<nrows; i++){
		lats[i] = 90.0 - i*180.0/(float)nrows;
	}
	for(i=0; i<ncols; i++){
		lons[i] = -180.0 + i*360/(float)ncols;
	}

	ret = nc_put_var_float(ncid, lat_varid, lats);
	assert(ret == NC_NOERR);
	ret = nc_put_var_float(ncid, lon_varid, lons);
	assert(ret == NC_NOERR);

	*_ncid = ncid;
	*_varid = varid;

	return 0;
}
