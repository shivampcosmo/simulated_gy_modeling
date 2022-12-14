#HDF5_C_INCLUDE=-I /usr/local/anaconda/anaconda3/include
#HDF5_C_LIBS=-L /usr/local/anaconda/anaconda3/lib -lpthread -lssl -lcrypto -lz -lm -lhdf5 -lhdf5_hl

HDF5_C_INCLUDE=-I /apps/python/3.6-conda5.2/include
HDF5_C_LIBS=-L /apps/python/3.6-conda5.2/lib -lpthread -lssl -lcrypto -lz -lm -lhdf5 -lhdf5_hl

CC = icc
CFLAGS = -Wall -O3 -ipo -qopenmp -std=c99
INCLDIRS = $(HDF5_C_INCLUDE) -I /usr/include -I ../
LFLAGS = $(HDF5_C_LIBS) -lgsl -lgslcblas # -lsvml
SOURCES = read_hdf5.c write_hdf5.c NFW_CDF.c hod.c compute_mocks.c
OBJECTS = $(SOURCES:.c=.o)
HEADERS = read_hdf5.h
EXEC = compute_mocks

.c.o:
	$(CC) $(CFLAGS) $(INCLDIRS) -c $<

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJECTS) $(LFLAGS)

$(OBJECTS): $(HEADERS) Makefile

clean:
	rm -f $(EXEC) *~ $(OBJECTS)


