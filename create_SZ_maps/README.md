# File details:
## create_LH_points.py:
This file will create the LH points for your parameters. You can change the params_all.ini file with providing a min and max values over which ther parameter needs to be varied. The format is that for each section you either give [min fiducial max] value for parameter to be varied, or you just give a [ fiducial ] value for it to be fixed. 
## clusters.py:
This is main code by Sigurd, opitmized to paint in Battaglia12 profile. I can get it to get a y-map from full sky websky catalog, cutting at halos above log(M) > 11.5 in ~1hr using 4 nodes of perlmutter. So generating O(100) simulated ymaps is not infeasible. 
## sim_websky_LH_ymaps.py:
This again is a modification of Sigurd's code, automated to readin the LH points and simply run a for-loop to get ymaps. It can also be paralellized, but since creating the main y-map is already parallelized and scales linearly with resources, one wouldn't gain much here I guess.


# Instructions on how to install and run on NERSC
## Installation on a fresh environment on perlmutter:
```
module load python
conda create --name myenv python=3.8
conda activate myenv

module load gsl
module load gcc

pip install healpy
git clone https://github.com/simonsobs/pixell.git
cd pixell
python setup.py build_ext -i

<!-- Replace this with your own pythonpath where pixell is downloaded -->
export PYTHONPATH=$PYTHONPATH:/global/cfs/cdirs/lsst/www/shivamp/pixell 

module load cray-fftw
module load cmake
pip install swig
export LDFLAGS+="-L$GSL_DIR/lib -L$FFTW_DIR"
export CPPFLAGS+="-I$GSL_DIR/include -I$FFTW_DIR/../include"
pip install pyccl

module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
MPICC="cc -shared" pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py

pip install camb
pip install smt
```

## Finally running the scripts:
### Get some interactive node to test it out:
<!-- Replace the account name with your nersc account, e.g. ACT or m1727 -->
```
salloc --nodes 2 --qos interactive --time 02:30:00 --constraint cpu --account=des
```

### Load everything:
```
module load python
conda activate myenv
module load gsl
module load gcc
export PYTHONPATH=$PYTHONPATH:/global/cfs/cdirs/lsst/www/shivamp/pixell
module load cray-fftw
module load cmake
export LDFLAGS+="-L$GSL_DIR/lib -L$FFTW_DIR"
export CPPFLAGS+="-I$GSL_DIR/include -I$FFTW_DIR/../include"
```
### Create the LH points:
```
python create_LH_points.py
```

### Run with MPI:
```
srun --nodes=2 --tasks-per-node=20 --cpu-bind=cores python sim_websky_LH_ymaps.py
```

