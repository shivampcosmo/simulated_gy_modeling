#include "read_hdf5.h"

#define M_PI 3.14159265358979323846
#define Mpc_to_cm 3.0856e24 /*Conversion factor from Mpc to cm */
#define Msun_to_g 1.989e33 /*Conversion factor from Msun to grams*/
#define G 6.672e-8 /*Universal Gravitational Constant in cgs units*/
#define Hubble 3.2407789e-18 /*Hubble's constant h/sec*/
#define rho_crit (3.0*pow(Hubble, 2.0) / (8.0 * M_PI * G)) * (pow(Mpc_to_cm, 3.0) / Msun_to_g) /*Cosmological Critical Density in Msun h^2 / Mpc^3 */
#define INDEX11(i,j) (i*11 + j)

/* These functions populate halos with HOD galaxies from an HDF5 halo catalog */
/* The outputs specify galaxy positions, and velocities (host mass and satellite/central identification have been dropped for speed and space considerations) */

hostDMH * find_galaxy_hosts(struct halo halos[], halo_metadata * env, double siglogM, double logMmin, double q_cen, double z_pivot, double slope, unsigned long int N_h, unsigned long int *Ncen, gsl_rng *r)
{
  /*This function uses the Zehavi 2011 prescription to find the halos that host central galaxies*/

  float * hosts = malloc(N_h * 11 * sizeof(float));
  int i;
  unsigned long j = 0;

  float f_logMmin_0 = (float)logMmin;
  float f_siglogM = (float)siglogM;
  float f_q_cen = (float)q_cen;

  for(i = 0; i < N_h; i++)
    {
      float logM = (float)log10(halos[i].mass);
      float env_rank = env[i].percentile;
      float prob;
	  
      /* compute modified logMmin for this halo's environment */
      float f_logMmin = f_logMmin_0 + slope * (halos[i].redshift - z_pivot) + f_q_cen * (env_rank - 0.5f);
      
      /*Mean central occupation or the probability of hosting a central*/

      prob = 0.5 * (1.0 + erf( (logM - f_logMmin) / f_siglogM) );
	  
      if(prob > gsl_rng_uniform(r))
	{
	  hosts[INDEX11(j,0)]  = halos[i].X;
	  hosts[INDEX11(j,1)]  = halos[i].Y;
	  hosts[INDEX11(j,2)]  = halos[i].Z;
	  hosts[INDEX11(j,3)]  = halos[i].vx;
	  hosts[INDEX11(j,4)]  = halos[i].vy;
	  hosts[INDEX11(j,5)]  = halos[i].vz;
	  hosts[INDEX11(j,6)]  = halos[i].mass;
	  hosts[INDEX11(j,7)]  = halos[i].redshift;   
	  hosts[INDEX11(j,8)]  = halos[i].conc;
	  hosts[INDEX11(j,9)]  = halos[i].R;
	  hosts[INDEX11(j,10)] = env[i].percentile;
	  j ++;
	}
    }
  
  *Ncen = j;

  hostDMH *host_coords = malloc(j*sizeof(hostDMH));

  for(i=0;i<j;i++)
    {
      host_coords[i].X              = hosts[INDEX11(i,0)];
      host_coords[i].Y              = hosts[INDEX11(i,1)];
      host_coords[i].Z              = hosts[INDEX11(i,2)];
      host_coords[i].vx             = hosts[INDEX11(i,3)];
      host_coords[i].vy             = hosts[INDEX11(i,4)];
      host_coords[i].vz             = hosts[INDEX11(i,5)];
      host_coords[i].mass           = hosts[INDEX11(i,6)];
      host_coords[i].redshift       = hosts[INDEX11(i,7)];
      host_coords[i].conc           = hosts[INDEX11(i,8)];
      host_coords[i].R              = hosts[INDEX11(i,9)];
      host_coords[i].env_percentile = hosts[INDEX11(i,10)];
    }

  return host_coords; 
}

int * find_satellites(struct halo halos[], double siglogM, double logMmin, double logM0, double logM1, double alpha, double q_sat, double z_pivot, double slope, unsigned long int N_h, unsigned long int *Nsat, gsl_rng *r)
{
  /*This function determines how many satellite galaxies each halo has using the same Zehavi 2011 prescription*/

  int i;
  unsigned long j =0;
  int * satellites = malloc(N_h*sizeof(int));

  for(i=0; i < N_h; i++)	 
    {
      double env_rank = (double)halos[i].env_percentile;
      double M0 = pow(10.0, logM0);
      double logM1_init = logM1;      
      double logM1 = logM1_init  + slope * (halos[i].redshift - z_pivot) + q_sat * (env_rank - 0.5f);
      double M1 = pow(10.0, logM1);	 
		      
      float logM = log10(halos[i].mass), M = halos[i].mass; //Come back to this once hdf5 is figured out
      double mean_cen = 1.0;
      double mean_sat;
      if(M < M0)
	{
	  mean_sat = 0.0; /* Enforcing the satellite cutoff */
	}
      else
	{
	  mean_sat = mean_cen * pow( ( ( M - M0 ) / M1 ), alpha );	  
	}
      unsigned int N_sat_this_halo = gsl_ran_poisson(r, mean_sat); /* Drawing from Poisson distribution centered at mean_sat to determine Nsat */
      satellites[i] =  N_sat_this_halo;
      j = j + N_sat_this_halo;
    }
  
  *Nsat = j;

  return satellites;
}

galaxy * pick_NFW_satellites(struct halo host, const int N_sat, double Acon, gsl_rng *r)
{
  /* This function determines the spatial distribution of the satellite galaxies */
  /* Galaxies are NFW-distributed using results from Correa et al. 2015 */
  /* We allow for a deviation from NFW in the galaxy distribution by a factor of r^delta_gamma */

  galaxy * coords = malloc(N_sat * sizeof(galaxy));

  double redshift = host.redshift;

  double alpha = 1.62774 - 0.2458*(1.0 + redshift) + 0.01716*pow(1.0 + redshift, 2.0);
  double beta = 1.66079 + 0.00359*(1.0 + redshift) - 1.6901*pow(1.0 + redshift, 0.00417);
  double gamma = -0.02049 + 0.0253*pow(1.0 + redshift ,-0.1044);
  	
  double logM = log10(host.mass), x0 = host.X, y0 = host.Y, z0 = host.Z, vx0 = host.vx, vy0 = host.vy, vz0 = host.vz;
  double exponent = alpha + beta*logM*(1.0 + gamma*pow(logM, 2.0)); /* Fit from Correa et al. 2015 */
  double cvir, Rvir;
  double M = host.mass, Grav = 4.301744 * pow(10, -9); /* Gravitational constant in (km/s)^2 (Mpc/Msun); G = 6.6743e-11 m^3 / ( kg s^2); (3.086e22)^-1 Mpc / m; 1.989e30 kg/Msun */
	
  cvir = Acon*sqrt(2.0)*pow(10.0, exponent); /* Approximate factor to rescale Rvir between crit, matter */
  Rvir = host.R;

  
  int j;

  /* pre-compute NFW profile for this satellite */
  float CDF[1000];
  size_t i;
  /* double prefac = (1.0 + cvir) / ( 2.0 + del_gamma - (1.0 + del_gamma)*gsl_sf_hyperg_2F1(1.0, 1.0, 3.0 + del_gamma, cvir/(1.0 + cvir)) ); */
  double prefac = 1.0 / ( log( 1.0 + cvir ) - (cvir / ( 1.0 + cvir )) ); /* Prefactor 1/A(c_vir) */
  float f_c_vir = (float)cvir;

#pragma simd
  for(i=0; i<1000; i++)
    {
      float x = (float)i / 1000.0;
      CDF[i] = prefac * ( log( 1.0 + x * f_c_vir ) - (x * f_c_vir / ( 1.0 + x*f_c_vir )) );
      /* CDF[i] = prefac * pow( x, (2.0 + del_gamma) ) * ( 2.0 + del_gamma - (1.0 + del_gamma)*gsl_sf_hyperg_2F1(1.0, 1.0, 3.0 + del_gamma, x*cvir/(x*cvir + 1.0)) ) / (1 + cvir*x); */
    }
  
  for(j=0; j<N_sat; j++)
    {
      double frac = NFW_CDF_sampler(&CDF[0], r);
      double R = Rvir * frac;
      double enc_mass = M * prefac * ( log( 1.0 + frac * f_c_vir ) - (frac * f_c_vir / ( 1.0 + frac*f_c_vir )) );
      double vel_disp = sqrt( (Grav * enc_mass) / (2.0*R*(1.0 / (1.0 + redshift)))); /*1D velocity disperion see eq. 14 Hikage, Takada, and Spergel 2018; factor of a = 1/1+z to convert from comoving to physical*/
      double phi = 2.0*M_PI*gsl_rng_uniform(r), costheta = 2.0*gsl_rng_uniform(r) - 1.0; /* Sphere point picking */
      double sintheta = sqrt(1.0 - costheta*costheta);
      double x = R*sintheta*cos(phi)+x0 , y = R*sintheta*sin(phi)+y0 , z = R*costheta+z0; /* Satellite Galaxy Coordinates */
      coords[j].X = x;
      coords[j].Y = y;
      coords[j].Z = z;
      coords[j].vx = vx0 + vel_disp * gsl_ran_gaussian(r, 1.0);
      coords[j].vy = vy0 + vel_disp * gsl_ran_gaussian(r, 1.0);
      coords[j].vz = vz0 + vel_disp * gsl_ran_gaussian(r, 1.0);
}

  return coords;
}

inline double wrap_periodic(double x, double Lbox)
{
  if((x < Lbox) && (x >= 0.)) {
    return x;
  } else if (x >= Lbox) {
    return (x-Lbox);
  } else if (x < 0.) {
    return (x+Lbox);
  }
}

void populate_hod(double siglogM, double logMmin, double logM0, double logM1, double alpha, double f_cen, double q_cen, double q_sat, double Acon, double z_pivot, double slope, unsigned long int seed, char *input_fname, char *output_fname, char *env_fname)
{
  herr_t status;
  size_t NumData,i;
  hostDMH *data;
  
  unsigned long Ncen;
  unsigned long Nsat;
  unsigned long fNcen = 0;
	
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed); /* Seeding random distribution */
	
  data = read_halo_hdf5(input_fname,"halos",&NumData);

  halo_metadata *env;
  size_t NumEnv;
  env = read_env_hdf5(env_fname,"halos",&NumEnv);
  	
  /* check that NumData == NumEnv */
  if(!(NumEnv == NumData)) exit(1);
	
  /* compute HOD parameters from number density, mass function, environment density */  
	
  hostDMH *cenhalos; //Central Coordinates
  cenhalos = find_galaxy_hosts(data, env, siglogM, logMmin, q_cen, z_pivot, slope, NumData, &Ncen, r);
  galaxy * cens = malloc(Ncen*sizeof(galaxy));
	
  for(i=0; i<Ncen; i++)
  {
    if(f_cen > gsl_rng_uniform(r))
    {
	double x = cenhalos[i].X, y = cenhalos[i].Y, z = cenhalos[i].Z;
	cens[fNcen].X = x;
      	cens[fNcen].Y = y;
      	cens[fNcen].Z = z;
	cens[fNcen].vx = cenhalos[i].vx;
      	cens[fNcen].vy = cenhalos[i].vy;
      	cens[fNcen].vz = cenhalos[i].vz;
        fNcen ++;
    }  
  }

  int *sats; //Number of Satellites for each halo
  sats = find_satellites(cenhalos, siglogM, logMmin, logM0, logM1, alpha, q_sat, z_pivot, slope, Ncen, &Nsat, r);
  galaxy * coords  = malloc(Nsat*sizeof(galaxy)); //Satellite Coordinates
  int j,k,l=0;
	
  for(j=0;j<Ncen;j++)
  {
    if(sats[j]>0){
      galaxy * halosats = malloc(sats[j] * sizeof(galaxy));
      halosats = pick_NFW_satellites(cenhalos[j], sats[j], Acon, r);
      for(k=0; k<sats[j]; k++)
	{ 
	  double x = halosats[k].X, y = halosats[k].Y, z = halosats[k].Z;
	  coords[l].X = x;
	  coords[l].Y = y;
	  coords[l].Z = z;
	  coords[l].vx = halosats[k].vx;
	  coords[l].vy = halosats[k].vy;
	  coords[l].vz = halosats[k].vz;
	  l++;
	}
      free(halosats);
    }
  }
  
  free(cenhalos);
  free(sats);
  	
  unsigned long int len = Nsat + fNcen;

  galaxy *HODgals = malloc(len*sizeof(galaxy));

  for(i=0; i<fNcen;i++)
    {
      HODgals[i].X = cens[i].X;
      HODgals[i].Y = cens[i].Y;
      HODgals[i].Z = cens[i].Z;
      HODgals[i].vx = cens[i].vx;	 
      HODgals[i].vy = cens[i].vy;	 
      HODgals[i].vz = cens[i].vz;
     }
  
  free(cens);
  	  
  for(i=0; i<Nsat; i++)
  {	  
      HODgals[i+fNcen].X = coords[i].X;
      HODgals[i+fNcen].Y = coords[i].Y;
      HODgals[i+fNcen].Z = coords[i].Z;
      HODgals[i+fNcen].vx = coords[i].vx;	 
      HODgals[i+fNcen].vy = coords[i].vy;	 
      HODgals[i+fNcen].vz = coords[i].vz;
  }
  free(coords);
	
  printf("Satellites Found. Writing to HDF5 file: %s\n", output_fname);
  status = write_gal_hdf5(output_fname, "particles", (size_t)len, HODgals);
  
  free(HODgals);
  
  gsl_rng_free(r);
}
