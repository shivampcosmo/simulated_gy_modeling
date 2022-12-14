#include "read_hdf5.h"

herr_t write_gal_hdf5(char filename[], char dataset_name[], size_t len, galaxy* data) {
  /* open HDF5 file*/
  hid_t memtype,filetype;
  hid_t file_id, dataset, space;
  hsize_t dims[1] = {len};
  herr_t status;

  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  memtype = H5Tcreate(H5T_COMPOUND, sizeof(galaxy));
  status = H5Tinsert(memtype, "x", HOFFSET(galaxy,X), H5T_NATIVE_FLOAT);
  status = H5Tinsert(memtype, "y", HOFFSET(galaxy,Y), H5T_NATIVE_FLOAT);
  status = H5Tinsert(memtype, "z", HOFFSET(galaxy,Z), H5T_NATIVE_FLOAT);  
  status = H5Tinsert(memtype, "vx", HOFFSET(galaxy,vx), H5T_NATIVE_FLOAT);
  status = H5Tinsert(memtype, "vy", HOFFSET(galaxy,vy), H5T_NATIVE_FLOAT);
  status = H5Tinsert(memtype, "vz", HOFFSET(galaxy,vz), H5T_NATIVE_FLOAT);

  size_t float_size_on_disk = H5Tget_size(H5T_IEEE_F32BE); // single precision
  size_t offset = 0;
  filetype = H5Tcreate(H5T_COMPOUND, 6*float_size_on_disk);
  status = H5Tinsert(filetype, "x", offset, H5T_IEEE_F32BE);
  offset += float_size_on_disk;
  status = H5Tinsert(filetype, "y", offset, H5T_IEEE_F32BE);
  offset += float_size_on_disk;
  status = H5Tinsert(filetype, "z", offset, H5T_IEEE_F32BE);
  offset += float_size_on_disk;
  status = H5Tinsert(filetype, "vx", offset, H5T_IEEE_F32BE);
  offset += float_size_on_disk;
  status = H5Tinsert(filetype, "vy", offset, H5T_IEEE_F32BE);
  offset += float_size_on_disk;
  status = H5Tinsert(filetype, "vz", offset, H5T_IEEE_F32BE);

  space = H5Screate_simple(1, dims, NULL);
  dataset = H5Dcreate(file_id, dataset_name, filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset,memtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

  H5Fclose(file_id);
  return status;
}
