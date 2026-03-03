#!/usr/bin/env python3



if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	from libstell.coils import COILSET
	from libstell.vmec import VMEC
	from libstell.bnorm import BNORM
	# from libstell.wall import WALL
	# from libstell.plot3D import PLOT3D
	# from libstell.libstell import FourierRep
	# import matplotlib.pyplot as pyplot
	# import numpy as np
	from stl import mesh
	parser = ArgumentParser(description= 
		'''Provides class for accessing coils files also serves as a
		   simple tool for assessing coils or coils files.''')
	parser.add_argument("-e", "--ext", dest="ext",
        help="Run all diagnostics for the given file extension", default = None)
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
	if args.ext:
		vmec_data = VMEC() 
		try:
			vmec_data.read_wout(args.ext)
		except:
			print(f'Could not file input file: wout_{args.vmec_ext}.nc or wout.{args.vmec_ext}')
			sys.exit(-1)
		coils.read_coils_file("coils."+args.ext)
		if args.mc:
			height,width = args.mc.split(',')
			coils_new = coils.multiToSingleFilament(int(height), int(width))
			coils = coils_new
		coils.plotcoils(lsave=args.lsaveplots)
		coils.print_coil_stats(lsave=args.lsaveplots)
		bnorm = BNORM()
		bnorm.read_bnorm_real("bnorm_real."+args.ext)
		bnorm.plot_bnorm_real_histogram(lsave=args.lsaveplots)
		bnorm.print_bnorm_real_stats()
		bnorm.read_bnorm_harm("bnorm_harm."+args.ext)
		bnorm.plot_bnorm_harm(type='sin', lsave=args.lsaveplots)
		bnorm.plot_bnorm_harm(type='cos', lsave=args.lsaveplots)
		coils.print_plasma_stats(bnorm.xreal, bnorm.yreal, bnorm.zzreal)

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
			bnorm.plot_bnorm_harm(type='cos', lsave=args.lsaveplots)
			