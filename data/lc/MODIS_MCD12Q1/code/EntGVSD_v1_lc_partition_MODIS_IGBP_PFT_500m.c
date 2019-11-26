
/************************************************************************
 * Author:  Qingsong Sun, Crystal Schaaf, 2016-05-01.                   *
 *                                                                      *
 * Ent Global Vegetation Structure Data set v1 processing step.         *
 *                                                                      *
 * Classifies the 500 m MODIS IGPB and PFT products into the 29 land    *
 * cover classes compatible with Ent PFTs, and outputs files in plain   *
 * binary format. Results are used as input files to produce coarser    *
 * grid maps with subgrid cover fractions of the above 29 classes.      *
 *                                                                      *
 * *********************************************************************/

#include <stdio.h>

#define N_COL (86400)
#define N_ROW (43200)

int main(int argc, char *argv[])
{
	if(argc != 4){
		printf("usage: %s IGBP.bin PFT.bin out.bin\n", argv[0]);
		return 1;
	}

	char *file_IGBP = argv[1];
	char *file_PFT  = argv[2];
	char *file_out  = argv[3];

	FILE *fp_IGBP = fopen(file_IGBP, "rb");
	if(fp_IGBP == NULL){
		printf("Open %s failed.\n", file_IGBP);
		return 1;
	}

	FILE *fp_PFT = fopen(file_PFT, "rb");
	if(fp_PFT == NULL){
		printf("Open %s failed.\n", file_PFT);
		fclose(fp_IGBP);
		return 1;
	}

	FILE *fp_out = fopen(file_out, "wb");
	if(fp_out == NULL){
		printf("Open %s failed.\n", file_out);
		fclose(fp_IGBP);
		fclose(fp_PFT);
		return 1;
	}


	unsigned char buf_IGBP[N_COL];
	unsigned char buf_PFT[N_COL];
	unsigned char buf_out[N_COL];

	int row, col;
	for(row=0; row<N_ROW; row++){
		printf("\r%3d%%", 100*row/N_ROW);
		fflush(stdout);

		if(N_COL != fread(buf_IGBP, 1, N_COL, fp_IGBP)){
			printf("Read %s failed. row=%d.\n", file_IGBP, row);
			fclose(fp_IGBP);
			fclose(fp_PFT);
			fclose(fp_out);
			return 1;
		}
		if(N_COL != fread(buf_PFT, 1, N_COL, fp_PFT)){
			printf("Read %s failed. row=%d.\n", file_PFT, row);
			fclose(fp_IGBP);
			fclose(fp_PFT);
			fclose(fp_out);
			return 1;
		}

		for(col=0; col<N_COL; col++){
			// init value
			//buf_out[col] = buf_IGBP[col];
			buf_out[col] = 28;
			
			// forest
			if (buf_IGBP[col] <= 4)
				buf_out[col] = buf_IGBP[col];

			//Mixed forest
			if ((buf_PFT[col] == 1 || buf_PFT[col] == 2) && buf_IGBP[col] == 5)
				buf_out[col] = 5;
			if (buf_PFT[col] >= 3 && buf_IGBP[col] == 5)
				buf_out[col] = 6;

			//closed shrublands
			if (buf_IGBP[col] == 6)
				buf_out[col] = 7;

			//open shrublands
			if (buf_IGBP[col] == 7)
				buf_out[col] = 8;

			// woody savannas
			if (buf_IGBP[col] == 8 && (buf_PFT[col] == 2 || buf_PFT[col] == 1))
				buf_out[col] = 9;

			if (buf_IGBP[col] == 8 && (buf_PFT[col] == 4 || buf_PFT[col] == 3))
				buf_out[col] = 10;

			if (buf_IGBP[col] == 8 && buf_PFT[col] == 5)
				buf_out[col] = 11;

			if (buf_IGBP[col] == 8 && buf_PFT[col] >= 6)
				buf_out[col] = 12;

			//savannas            
			if (buf_IGBP[col] == 9 && buf_PFT[col] <= 5 && buf_PFT[col] > 0)
				buf_out[col] = 13;

			if (buf_IGBP[col] == 9 && buf_PFT[col] >= 6)
				buf_out[col] = 14;

			// Grassland

			if (buf_IGBP[col] == 10)
				buf_out[col] = 15;

			// permanent wetlands
			if (buf_IGBP[col] == 11 && buf_PFT[col] <= 4 && buf_PFT[col] > 0)
				buf_out[col] = 16;

			//if(buf_IGBP[col]==11&&(buf_PFT[col]==3 ||buf_PFT[col]==4))
			//buf_out[col]=17;

			if (buf_IGBP[col] == 11 && buf_PFT[col] >= 5)
				buf_out[col] = 17;

			//if(buf_IGBP[col]==11&&buf_PFT[col]>=6)
			//buf_out[col]=19;

			// croplands

			if (buf_IGBP[col] == 12 && buf_PFT[col] <= 4 && buf_PFT[col] > 0)
				buf_out[col] = 19;

			if (buf_IGBP[col] == 12 && (buf_PFT[col] == 5 || buf_PFT[col] == 8))
				buf_out[col] = 19;

			if (buf_IGBP[col] == 12 && (buf_PFT[col] == 6 || buf_PFT[col] == 7))
				buf_out[col] = 18;

			if (buf_IGBP[col] == 12 && buf_PFT[col] > 8)
				buf_out[col] = 19;

			// urban_built-up
			if (buf_IGBP[col] == 13)
				buf_out[col] = 20;

			//cropland/natural vegetation

			if (buf_IGBP[col] == 14 && buf_PFT[col] <= 2 && buf_PFT[col] > 0)
				buf_out[col] = 21;

			if (buf_IGBP[col] == 14 && (buf_PFT[col] == 3 || buf_PFT[col] == 4))
				buf_out[col] = 22;

			if (buf_IGBP[col] == 14 && buf_PFT[col] == 5)
				buf_out[col] = 23;

			if (buf_IGBP[col] == 14 && buf_PFT[col] == 6)
				buf_out[col] = 24;

			if (buf_IGBP[col] == 14 && buf_PFT[col] == 7)
				buf_out[col] = 25;

			if (buf_IGBP[col] == 14 && buf_PFT[col] >= 8)
				buf_out[col] = 26;

			// permanent snow and ice

			if (buf_IGBP[col] == 15)
				buf_out[col] = 27;

			// barren or sparely vegetated

			if (buf_IGBP[col] == 254){
				buf_out[col] = 28;
			}	
		}
	
		// write output
		if(N_COL != fwrite(buf_out, 1, N_COL, fp_out)){
			printf("write %s failed. row=%d.\n", file_out, row);
			fclose(fp_IGBP);
			fclose(fp_PFT);
			fclose(fp_out);
			return 1;
		}
	}
	
	fclose(fp_IGBP);
	fclose(fp_PFT);
	fclose(fp_out);
	return 0;
}
