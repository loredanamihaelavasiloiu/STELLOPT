#!/usr/bin/env python3

if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.coils import COILSET
	from libstell.vmec import VMEC
	from libstell.bnorm import BNORM
	from libstell.fieldlines import FIELDLINES
	from libstell.stellopt import STELLOPT
	# from libstell.wall import WALL
	# from libstell.plot3D import PLOT3D
	# from libstell.libstell import FourierRep
	import matplotlib.pyplot as pyplot
	import numpy as np
	from stl import mesh
	parser = ArgumentParser(description= 
		'''Provides class for accessing coils files also serves as a
		   simple tool for assessing coils or coils files.''')
	parser.add_argument("-e", "--ext", dest="ext",
        help="STELLOPT files extension", default = None)
	parser.add_argument("-ll", "--loglevel", dest="loglevel",
		help="Logging level, where higher means more output. Negative loglevel means exclusively that level, positive includes all lower levels." \
		" Options: -3, -2, -1, 0, 1, 2, 3. Default: 0", default = 0, type=int)
	parser.add_argument("-mc", dest="mc",
        help="Signals that the input coil file is a multi-filament coil and should be averaged to a single filament coil", default = None)
	
	parser.add_argument("-c", "--coil", dest="coils_file",
		help="Coils file for input", default = None)
	
	parser.add_argument("-v", "--vmec", dest="vmec_ext",
		help="VMEC file extension.", default = None)
	
	parser.add_argument("-br", "--bnorm_real", dest="bnormreal_file",
        help="bnorm_real file for input", default = None)	
	
	parser.add_argument("-bh", "--bnorm_harm", dest="bnormharm_file",
        help="bnorm_harm file for input", default = None)

	parser.add_argument("-f", "--fieldlines", dest="fieldlines_ext",
		help="Field lines file for input. Required for plotting field lines.", default = None)

	parser.add_argument("-pfl", "--plot_fieldlines", dest="lplotfieldlines", action='store_true',
		help="Plot the field lines output. Requires FIELDLINES output", default = False)

	parser.add_argument("-cp", "--coils_plot", dest="lplotcoil", action='store_true',
		help="Plot the coils file.", default = False)
	
	parser.add_argument("-cs", "--coils_stats", dest="lcoil_stats", action='store_true',
        help="Print some statistics of the coils and coil-to-coil clearances.", default = False)
	
	parser.add_argument("-cp_stats", "--coils_plasma_stats", dest="lcoil_plasmastats", action='store_true',
        help="Print some statistics of the coil-to-plasma clearance. Requires bnorm_real and coils file.", default = False)
	
	parser.add_argument("--bnorm_real_histogram", dest="lbnormreal_histogram", action='store_true',
        help="Plot a histogram of the bnorm values.", default = False)
	
	parser.add_argument("--bnorm_real_stats", dest="lbnormreal_stats", action='store_true',
        help="Print some statistics of the bnorm values.", default = False)
	
	parser.add_argument("--bnorm_harm_plot", dest="lbnormharm_plot", action='store_true',
        help="Plot the bnorm_harm values.", default = False)
	
	parser.add_argument("--saveplots", dest="lsaveplots", action='store_true',
        help="Save the generated plots.", default = False)
	
	args = parser.parse_args()
	coils = COILSET()
	if args.lsaveplots:
		import os
		# Check if the directory with the name of the current ext exists, if not create it
		if not os.path.exists(args.ext):
			os.makedirs(args.ext)
			print(f'Created directory: {args.ext}')
		else:
			print(f'Directory already exists: {args.ext}')
	if args.ext:
		stellopt = STELLOPT()
		vmec_data = VMEC() 
		bnorm = BNORM()
		try:
			vmec_data.read_wout(args.ext)
		except:
			print(f'Could not file input file: wout_{args.ext}.nc or wout.{args.ext}')
			sys.exit(-1)
		coils.read_coils_file(args.ext)
		if args.mc:
			height,width = args.mc.split(',')
			coils_new = coils.multiToSingleFilament(int(height), int(width))
			coils = coils_new
		bnorm.read_bnorm('bnorm.'+args.ext)
		bnorm.read_bnorm_real(args.ext)
		bnorm.read_bnorm_harm(args.ext)
		if args.loglevel > 0 or args.loglevel == -1:
			coils.plotcoils(lsave=args.lsaveplots)
			coils.print_coil_stats(lsave=args.lsaveplots)
			bnorm.plot_bnorm_real_histogram(lsave=args.lsaveplots)
			bnorm.print_bnorm_real_stats()
			bnorm.plot_bnorm_harm(type='sin', lsave=args.lsaveplots)
#			bnorm.plot_bnorm_harm(type='cos', lsave=args.lsaveplots)
			coils.print_plasma_stats(bnorm.xreal, bnorm.yreal, bnorm.zzreal, lsave=args.lsaveplots)
			from pathlib import Path
			if not Path(f'fieldlines_{args.ext}.h5').is_file():
				print(f'Field lines file not found: fieldlines_{args.ext}.h5. Run FIELDLINES to generate this file before plotting field lines.')
			else:
				field_data = FIELDLINES()
				field_data.read_fieldlines('fieldlines_'+args.ext+'.h5')
				phi0 = 0
				field_data.plot_poincare(0, lsave=args.lsaveplots)
				phi1 = field_data.PHI_lines[0,int(np.round(field_data.npoinc/4))]
				field_data.plot_poincare(phi1, lsave=args.lsaveplots)
				phi2 = field_data.PHI_lines[0,int(np.round(field_data.npoinc/2))]
				field_data.plot_poincare(phi2, lsave=args.lsaveplots)
		if args.loglevel > 1 or args.loglevel == -2:
			# BNORM values on the plasma surface
			bnorm.plotBrealsurf(lsave=args.lsaveplots)
			# RMS(bnorm total) vs. iteration number
			# Max(bnorm_total) vs iteration number
			# "dominant" Fourier mode of bnorm_harm vs iteration number

			# Read coil curvature file
			stellopt.read_stellopt_coil_curvature('coil_curvature.'+args.ext)
			# Coil curvature diagnostics (profile, max, rms, histogram)
			stellopt.print_coil_curvature_stats()
			stellopt.plot_stellopt_coil_curvature(lsave=args.lsaveplots, save_folder=args.ext)
			# Coil torsion diagnostics (profile, max, rms, histogram)
			stellopt.print_coil_torsion_stats()
			stellopt.plot_stellopt_coil_torsion(lsave=args.lsaveplots, save_folder=args.ext)
			# Total objective vs iteration number
			# Individual objective terms vs iteration number
			# Optimization variables vs iteration number
			# Coil magnetic energy
			pass

		if args.loglevel > 2 or args.loglevel == -3:
			# BNORM plasma vs coil vs total comparison
			# Relative BNORM error
			# Axis-field diagnostic (B.T along magnetic axis)
			# Coils + plasma overlay plot
			# Arclength discretization quality check
			# --- Spline knot 3D plot (geometry) with method for visualizing geometric regions bound by x_coil_kts_min and x_coil_kts_max
			pass


	if args.vmec_ext:
		vmec_data = VMEC()
		try:
			vmec_data.read_wout(args.vmec_ext)
		except:
			print(f'Could not file input file: wout_{args.vmec_ext}.nc or wout.{args.vmec_ext}')
			sys.exit(-1)

	if args.coils_file: 
		coils.read_coils_file(args.coils_file)
		if args.lplotcoil: coils.plotcoils(lsave=args.lsaveplots)
		if args.lcoil_stats: coils.print_coil_stats(lsave=args.lsaveplots)

	if args.bnormreal_file:
		bnorm = BNORM()
		bnorm.read_bnorm_real(args.bnormreal_file)
		if args.lbnormreal_histogram:
			bnorm.plot_bnorm_real_histogram(lsave=args.lsaveplots)
		if args.lbnormreal_stats:
			bnorm.print_bnorm_real_stats()
		
		if args.coils_file and args.lcoil_plasmastats:
			coils.print_plasma_stats(bnorm.xreal, bnorm.yreal, bnorm.zzreal, lsave=args.lsaveplots)
        

	if args.bnormharm_file:
		bnorm = BNORM()
		bnorm.read_bnorm_harm(args.bnormharm_file)
		if args.lbnormharm_plot:
			bnorm.plot_bnorm_harm(type='sin', lsave=args.lsaveplots)
#			bnorm.plot_bnorm_harm(type='cos', lsave=args.lsaveplots)

	if args.fieldlines_ext:
		field_data = FIELDLINES()
		field_data.read_fieldlines('fieldlines_'+args.fieldlines_ext+'.h5')
		if args.lplotfieldlines:
			phi0 = 0
			field_data.plot_poincare(0, lsave=args.lsaveplots)
			phi1 = field_data.PHI_lines[0,int(np.round(field_data.npoinc/4))]
			field_data.plot_poincare(phi1, lsave=args.lsaveplots)
			phi2 = field_data.PHI_lines[0,int(np.round(field_data.npoinc/2))]
			field_data.plot_poincare(phi2, lsave=args.lsaveplots)
			