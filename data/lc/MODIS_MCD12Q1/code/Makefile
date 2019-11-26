#
TARGET1 = EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_500m.exe
TARGET2 = EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_subgrid_1km.exe
TARGET3 = merge.exe

# Files
OBJ1 = EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_500m.o
OBJ2 = EntGVSD_v1_lc_partition_MODIS_IGBP_PFT_subgrid_1km.o
OBJ3 = merge.o

##########################################3
ADD_CFLAGS= -O3 -DLYNX  -ffloat-store

INC = -I. -I$(NC_INC) -I$(HDFINC) -I$(HDFEOS_INC)

LIB = -lz -lm -L$(HDFLIB) -ldf -lmfhdf -L$(NC_LIB) -lnetcdf  -L$(HDFEOS_LIB) -lhdfeos -lGctp

ALL : $(TARGET1) $(TARGET2) $(TARGET3)

# make
$(TARGET1) : $(OBJ1)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(OBJ1) $(LIB) -o $(TARGET1)

$(TARGET2) : $(OBJ2)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(OBJ2) $(LIB) -o $(TARGET2)

$(TARGET3) : $(OBJ3)
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(OBJ3) $(LIB) -o $(TARGET3)


.c.o: 
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INC) -c $< -o $@

#delete object files:
clean:
	rm *.o
	rm *.exe

