# This is a fork from Sigurd's code. The goal here is simply to automate it
# to run on a latin hypercube grid. One can vary the pressure profile parameters.
# Note that for the profile shape, only beta parameter can be varied. alpha and gamma are fixed.

# From Sigurd:
# The goal of this program is to read in the WebSky cluster catalog
# and produce an enmap for the given geometry, frequency and beam.
# We want to reproduce what Nemo does, but with my cluster profile
# evaluation code

from pixell import utils, bunch
import numpy as np, pyccl, time
import configparser
# config = configparser.ConfigParser()

fname_fid_params = 'params_fid.ini'
do_run_fid_only = False
nparams = 3
rsv = 0
# config.read(fname_fid_params)
if not do_run_fid_only:
	LH_points_fname = 'sample_chain_vary_params_np_' + str(nparams) + '_rs' + str(rsv) + '.txt'
	LH_points = np.loadtxt(LH_points_fname)
	with open(LH_points_fname) as f:
		LH_header = f.readline()
	LH_header_split = LH_header[1:].split()
	params_split = LH_header_split[1:]
	nLH = len(LH_points)
else:
	nLH = 1
	print("Running fiducial only")

for jb in range(nLH):
	print('DOING JOB: ',jb)
	config_fid_params = configparser.ConfigParser()	
	config_fid_params.read(fname_fid_params)

	if not do_run_fid_only:
		LH_jb = LH_points[jb][1:]
		for jp in range(len(params_split)):
			sec = params_split[jp].split('--')[0]
			var = params_split[jp].split('--')[1]
			config_fid_params[sec][var] = str(LH_jb[jp])
	
	pressure_all_jb = config_fid_params['pressure_params']
	
	P0_A, P0_alm, P0_alz = float(pressure_all_jb['p0_a_m']), float(pressure_all_jb['p0_alpha_m']), float(pressure_all_jb['p0_alpha_z'])
	xc_A, xc_alm, xc_alz = float(pressure_all_jb['xc_a_m']), float(pressure_all_jb['xc_alpha_m']), float(pressure_all_jb['xc_alpha_z'])
	beta_A, beta_alm, beta_alz = float(pressure_all_jb['beta_a_m']), float(pressure_all_jb['beta_alpha_m']), float(pressure_all_jb['beta_alpha_z'])
	bpl_alm = float(pressure_all_jb['bpl_alm'])
	Mbreak = 1e14 # Msun/h
	param_updates_all = bunch.Bunch(P0_A=P0_A, P0_alm=P0_alm, P0_alz=P0_alz, 
	xc_A=xc_A, xc_alm=xc_alm, xc_alz=xc_alz, beta_A=beta_A, 
	beta_alm=beta_alm, beta_alz=beta_alz, bpl_alm=bpl_alm, Mbreak=Mbreak)
	broken_pl = False
	import os.path
	if not do_run_fid_only:
		save_ldir = 'output_ymaps/abacus_nLHtot_' + str(nLH) + '_np_' + str(nparams) + '_rs_' + str(rsv) 
		#open folder if does not exist
		if not os.path.exists(save_ldir):
			os.makedirs(save_ldir)
		else:
			pass
		save_fname =  save_ldir + '/sim_abacus_nLHtot_' + str(nLH) + '_np_' + str(nparams) + '_rs_' + str(rsv) +  '_LHid_' + str(jb) + '.fits'
	else:
		save_fname = 'output_ymaps/sim_abacus_fid_B12_test.fits'
	# check if file does not exist
	if os.path.isfile(save_fname):
		pass
	else:
		import argparse
		parser = argparse.ArgumentParser()
		parser.add_argument("-halos", "--halos", type=str, default='data_inp/abacus_halo_cat_Mmin_1e11.fits')
		parser.add_argument("-geometry", "--geometry", type=str, default='data_inp/abacus_geometry.fits')
		parser.add_argument("-ofile", "--ofile", type=str, default=save_fname,   help="Output map")
		parser.add_argument("-f", "--freq",  type=float, default=98, help="Frequency to simulate at")
		parser.add_argument("-b", "--beam",  type=str,   default="1.6", help="The beam to use. Either a number, which is interpreted as a FWHM in arcmin, or a file name, which should be to a beam transform file with the format [l,b(l)]")
		parser.add_argument("-n", "--nhalo", type=int,   default=0, help="Number of halos to use. 0 for unlimited")
		parser.add_argument("-lg10Mmin", "--lg10Mmin", type=float,   default=11.5, help="Minimum halo mass to use. 11.5 should suffice for cross-correlation purposes")
		parser.add_argument("-B", "--bsize", type=int,   default=1000)
		parser.add_argument("-V", "--vmin",  type=float, default=0.001)
		args = parser.parse_args()

		from astropy.io import fits
		from pixell import enmap, utils, bunch, pointsrcs
		import clusters
		from mpi4py import MPI as mpi

		dtype   = np.float32 # only float32 supported by fast srcsim
		comm    = mpi.COMM_WORLD
		bsize   = args.bsize
		margin  = 100
		beta_range = [-14,-3]
		# This is the profile cutoff in µK. It defaults to 1e-3, i.e. 1 nK. This
		# should be much lower than necessary, but I set it this low because the
		# profile building step is the bottleneck for most objects anyway. If that
		# step could be sped up, then this program could be made much faster by increasing
		# this number. But for now it's practically free to keep it this low.
		cosmology   = pyccl.Cosmology(Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159, n_s=0.9667, transfer_function="boltzmann_camb")
		vmin    = args.vmin
		lg10Mmin = args.lg10Mmin
		catalog_data = fits.open(args.halos)[1].data
		Mhalo_all = catalog_data['M200c']
		indsel = np.where((Mhalo_all > 10**lg10Mmin))[0]
		nhalo   = len(indsel)
		z = catalog_data['redshift']
		ra, dec = utils.rect2ang(np.array([catalog_data['X'], catalog_data['Y'], catalog_data['Z']]))
		# ra = catalog_data['ra']
		# dec = catalog_data['dec']
		m200 = np.array(catalog_data['M200c'], dtype=np.float64)
		if args.lg10Mmin is not None:
			Mmin = (10**lg10Mmin)
			indsel = np.where((m200 > Mmin))[0]
			z_all, ra_all, dec_all, m200_all = z[indsel], ra[indsel], dec[indsel], (m200[indsel])*float(utils.M_sun/cosmology["h"])
		# import pdb; pdb.set_trace()
		# if args.nhalo:
		# 	nhalo = min(nhalo, args.nhalo)

		
		try:
			shape, wcs  = enmap.read_map_geometry(args.geometry)
		except:
			# shape,wcs = enmap.fullsky_geometry(res=1.0 * utils.arcmin)
			dec_from, dec_to = np.deg2rad([0,90.])
			ra_from, ra_to = np.deg2rad([0, 90])
			box = [[dec_from,ra_from],[dec_to,ra_to]]
			shape,wcs = enmap.geometry(pos=box,res=1.0 * utils.arcmin)
			# shape,wcs = enmap.fullsky_geometry(res=1.0 * utils.arcmin)		
		freq        = args.freq*1e9
		omap        = enmap.zeros(shape[-2:], wcs, dtype)
		rht         = utils.RadialFourierTransform()
		# Read the beam from one of the two formats
		try:
			sigma = float(args.beam)*utils.fwhm*utils.arcmin
			lbeam = np.exp(-0.5*rht.l**2*sigma**2)
		except ValueError:
			l, bl = np.loadtxt(args.beam, usecols=(0,1), ndmin=2).T
			bl   /= np.max(bl)
			lbeam = np.interp(rht.l, l, bl)
		prof_builder= clusters.ProfileBattagliaFast(cosmology=cosmology, beta_range=beta_range, params_update=param_updates_all, 
				broken_pl=broken_pl)
		mass_interp = clusters.MdeltaTranslator(cosmology)

		# We use this to decide if it's worth it to prune
		# objects outide our patch or now. This pruning takes
		# some extra calculations which aren't necessary if we're
		# fullsky or close to it
		fullsky = enmap.area(shape, wcs)/(4*np.pi) > 0.8

		# Loop over halos
		nblock  = (nhalo+bsize-1)//bsize
		tget    = 0
		tprof   = 0
		tpaint  = 0
		ntot    = 0
		for bi in range(comm.rank, nblock, comm.size):
			i1    = bi*bsize
			i2    = min((bi+1)*bsize, nhalo)
			t1    = time.time()
			# data = 
			z_sel, ra_sel, dec_sel, m200_sel = z_all[i1:i2], ra_all[i1:i2], dec_all[i1:i2], m200_all[i1:i2]
			print(i1, i2)
			# data  = clusters.websky_pkcs_read(args.halos, num=i2-i1, offset=i1)
			# Prune the ones outside our area
			# if not fullsky:
			# 	pixs  = enmap.sky2pix(shape, wcs, utils.rect2ang(data.T[:3])[::-1])
			# 	good  = np.all((pixs >= -margin) & (pixs < np.array(shape[-2:])[:,None]+margin), 0)
			# 	data  = data[good]
			ngood = len(z_sel)
			if ngood == 0: continue
			# Compute physical quantities
			# cat    = clusters.websky_decode(data, cosmology, mass_interp, args.lg10Mmin); del data
			cat = bunch.Bunch(z=z_sel, ra=ra_sel, dec=dec_sel, m200=m200_sel)
			del z_sel, ra_sel, dec_sel, m200_sel
			t2     = time.time(); tget = t2-t1
			# Evaluate the y profile
			rprofs  = prof_builder.y(cat.m200[:,None], cat.z[:,None], rht.r)
			# convolve with beam
			lprofs  = rht.real2harm(rprofs)
			lprofs  *= lbeam
			rprofs  = rht.harm2real(lprofs)
			r, rprofs = rht.unpad(rht.r, rprofs)
			# and factor out peak value
			yamps   = rprofs[:,0].copy()
			rprofs /= yamps[:,None]
			# Prepare for painting
			amps   = (yamps * utils.tsz_spectrum(freq) / utils.dplanck(freq) * 1e6).astype(dtype)
			poss   = np.array([cat.dec,cat.ra]).astype(dtype)
			profiles = [np.array([r,prof]).astype(dtype) for prof in rprofs]; del rprofs
			prof_ids = np.arange(len(profiles)).astype(np.int32)
			# And paint
			ntot += ngood
			t3 = time.time(); tprof  = t3-t2
			pointsrcs.sim_objects(shape, wcs, poss, amps, profiles, prof_ids=prof_ids, omap=omap, vmin=vmin)
			t4 = time.time(); tpaint = t4-t3
			# Print our status
			print("%3d %3d %4d/%d ndone %6.1fk get %6.3f ms prof %6.3f ms draw %6.3f tot %6.3f each maxamp %7.2f" % (
				jb, comm.rank, bi+1, nblock, ntot/1e3, tget/ngood*1e3, tprof/ngood*1e3, tpaint/ngood*1e3,
				(tget+tprof+tpaint)/ngood*1e3, np.max(np.abs(amps))))
			# import pdb; pdb.set_trace()

		print("%4d Reducing" % comm.rank)

		if comm.size > 1:
			comm.Barrier()
			if comm.rank == 0:
				omap_full = np.zeros_like(omap)
				comm.Reduce(omap, omap_full)
				omap = omap_full
			else:
				comm.Reduce(omap, None)
				del omap
			comm.Barrier()
		if comm.rank == 0:
			enmap.write_map(args.ofile, omap)
			print("Done")
