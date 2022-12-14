#include "read_hdf5.h"

int main(int argc, char *argv[])
{
  if(argc != 16) {
    printf("%d usage: ./compute_mocks siglogM logMmin logM0 logM1 alpha f_cen q_cen q_sat Acon z_pivot slope seed [halo catalog file] [galaxy mock file] [halo environment file]\n",argc);
    return -1;
  }

  double siglogM = strtod(argv[1], NULL);
  double logMmin = strtod(argv[2], NULL);
  double logM0 = strtod(argv[3], NULL);
  double logM1 = strtod(argv[4], NULL);
  double alpha = strtod(argv[5], NULL);
  double f_cen = strtod(argv[6], NULL);
  double q_cen = strtod(argv[7], NULL);
  double q_sat = strtod(argv[8], NULL);
  double Acon = strtod(argv[9], NULL);
  double z_pivot = strtod(argv[10], NULL);
  double slope = strtod(argv[11], NULL);		
 	
  int seed = atoi(argv[12]);

  char *halo_file, *output_file, *env_file;
  size_t halo_file_ssize, output_file_ssize, env_file_ssize;

  halo_file_ssize = sizeof(char)*(strlen(argv[13]) + 1);
  output_file_ssize = sizeof(char)*(strlen(argv[14]) +1 );
  env_file_ssize = sizeof(char)*(strlen(argv[15]) +1 );

  halo_file = malloc(halo_file_ssize);
  output_file = malloc(output_file_ssize);
  env_file = malloc(env_file_ssize);
  snprintf(halo_file, halo_file_ssize, "%s", argv[13]);
  snprintf(output_file, output_file_ssize, "%s", argv[14]);
  snprintf(env_file, env_file_ssize, "%s", argv[15]);

  fprintf(stderr,"Computing HOD from %s\n", halo_file);
  fprintf(stderr,"Reading environment density from %s\n", env_file);
  fprintf(stderr,"Saving to output file %s\n", output_file);

  populate_hod(siglogM, logMmin, logM0, logM1, alpha, f_cen, q_cen, q_sat, Acon, z_pivot, slope,\
	       seed, halo_file, output_file, env_file);

  return 0;
}
