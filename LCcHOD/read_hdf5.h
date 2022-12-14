#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <time.h>

#include "hdf5.h"
#include "hdf5_hl.h"

typedef struct halo
{
  float mass;
  float X;
  float Y;
  float Z;
  float vx;
  float vy;
  float vz;
  float conc;
  float R;
  float redshift;
  float env_percentile;
} hostDMH;

typedef struct HODgal
{
  float X;
  float Y;
  float Z;
  float vx;
  float vy;
  float vz;
} galaxy;

typedef struct _halo_metadata
{
  float mass;
  float density;
  float percentile;
} halo_metadata;

void populate_hod(double siglogM, double logMmin, double logM0, double logM1, double alpha, double f_cen, double q_cen, double q_sat, double Acon, double z_pivot, double slope, unsigned long int seed, char *input_fname, char *output_fname, char *env_fname);
double NFW_CDF_sampler(float * restrict CDF, gsl_rng *r);

void* read_halo_hdf5(char infile[],char dataset_name[],size_t *len);
herr_t write_gal_hdf5(char filename[], char dataset_name[], size_t len, galaxy* data);
void* read_env_hdf5(char filename[], char dataset_name[], size_t *len);
