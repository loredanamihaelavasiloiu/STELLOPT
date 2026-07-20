##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling BNORM
data.
"""

# Libraries
from libstell.libstell import LIBSTELL, FourierRep

# Constants

# BNORM Class
class BNORM(FourierRep):
	"""Class for working with BNORM equilibria

	"""
	def __init__(self):
		super().__init__()
		self.libStell = LIBSTELL()

	def read_bnorm(self,filename):
		"""Reads a BNORM file

		This routine reads and initilizes the BNORM class
		with variable information from a BNORM file.

		Parameters
		----------
		file : str
			Path to wout file.
		"""
		import numpy as np
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		self.mnmax = len(lines)
		self.xm = np.zeros((self.mnmax,1))
		self.xn = np.zeros((self.mnmax,1))
		self.bnmnc = np.zeros((1,self.mnmax))
		self.bnmns = np.zeros((1,self.mnmax))
		mn = 0
		for line in lines:
			(txt1,txt2,txt3) = line.split()
			self.xm[mn] = int(txt1)
			self.xn[mn] = int(txt2)
			self.bnmns[0,mn] = float(txt3)
			mn = mn + 1

	def read_bnorm_real(self,ext):
		"""Reads a BNORM_REAL file

		This routine reads and initilizes the BNORM class
		with variable information from a BNORM_REAL file.

		Parameters
		----------
		filename : str
			Path to bnorm_real file.
		"""
		import numpy as np
		self.ext = ext
		f = open("bnorm_real."+ext,'r')
		lines = f.readlines()
		f.close()
		self.nuv=int(lines[0])
		u             = np.zeros(self.nuv,dtype=np.int64)
		v             = np.zeros(self.nuv,dtype=np.int64)
		theta         = np.zeros(self.nuv)
		zeta          = np.zeros(self.nuv)
		phi           = np.zeros(self.nuv)
		rreal         = np.zeros(self.nuv)
		zreal         = np.zeros(self.nuv)
		Nx            = np.zeros(self.nuv)
		Ny            = np.zeros(self.nuv)
		Nz            = np.zeros(self.nuv)
		bnreal        = np.zeros(self.nuv)
		bcreal        = np.zeros(self.nuv)
		bnormal_total = np.zeros(self.nuv)
		for j in range(self.nuv):
			txt                   = lines[j+1].split()
			u[j]             = int(txt[1])
			v[j]             = int(txt[2])
			theta[j]         = float(txt[3])
			zeta[j]          = float(txt[4])
			phi[j]           = float(txt[5])
			rreal[j]         = float(txt[6])
			zreal[j]         = float(txt[7])
			Nx[j]            = float(txt[8])
			Ny[j]            = float(txt[9])
			Nz[j]            = float(txt[10])
			bnreal[j]        = float(txt[11])
			bcreal[j]        = float(txt[12])
			bnormal_total[j] = float(txt[13])

		# Field periodicity of the stellarator
		ratio = phi/zeta
		# Remove Nans from ratio
		ratio = ratio[~np.isnan(ratio)]
		Nfp = int(2/np.mean(ratio))

		self.xreal = np.zeros(self.nuv*Nfp)
		self.yreal = np.zeros(self.nuv*Nfp)
		self.zzreal = np.zeros(self.nuv*Nfp)
		self.xreal[0:self.nuv] = rreal * np.cos(phi)
		self.yreal[0:self.nuv] = rreal * np.sin(phi)
		self.zzreal[0:self.nuv] = zreal
		self.xreal[self.nuv:2*self.nuv] = rreal * np.sin(phi)
		self.yreal[self.nuv:2*self.nuv] = rreal * np.cos(phi)
		self.zzreal[self.nuv:2*self.nuv] = -zreal
		self.xreal[2*self.nuv:3*self.nuv] = -rreal * np.sin(phi)
		self.yreal[2*self.nuv:3*self.nuv] = rreal * np.cos(phi)
		self.zzreal[2*self.nuv:3*self.nuv] = zreal
		self.xreal[3*self.nuv:4*self.nuv] = -rreal * np.cos(phi)
		self.yreal[3*self.nuv:4*self.nuv] = rreal * np.sin(phi)
		self.zzreal[3*self.nuv:4*self.nuv] = -zreal
		self.xreal[4*self.nuv:5*self.nuv] = rreal * np.cos(phi)
		self.yreal[4*self.nuv:5*self.nuv] = -rreal * np.sin(phi)
		self.zzreal[4*self.nuv:5*self.nuv] = -zreal
		self.xreal[5*self.nuv:6*self.nuv] = rreal * np.sin(phi)
		self.yreal[5*self.nuv:6*self.nuv] = -rreal * np.cos(phi)
		self.zzreal[5*self.nuv:6*self.nuv] = zreal
		self.xreal[6*self.nuv:7*self.nuv] = -rreal * np.sin(phi)
		self.yreal[6*self.nuv:7*self.nuv] = -rreal * np.cos(phi)
		self.zzreal[6*self.nuv:7*self.nuv] = -zreal
		self.xreal[7*self.nuv:8*self.nuv] = -rreal * np.cos(phi)
		self.yreal[7*self.nuv:8*self.nuv] = -rreal * np.sin(phi)
		self.zzreal[7*self.nuv:8*self.nuv] = zreal
		# for i in range(Nfp):
			# self.xreal[i*self.nuv:(i+1)*self.nuv] = rreal * np.cos(phi)
			# self.yreal[i*self.nuv:(i+1)*self.nuv] = rreal * np.sin(phi)
			# self.zzreal[i*self.nuv:(i+1)*self.nuv] = zreal
		self.full_bnormal_total = np.tile(bnormal_total, Nfp)

		nu = max(u)
		nv = max(v)
		self.u             = u.reshape(nu,nv)
		self.v             = v.reshape(nu,nv)
		self.theta          = theta.reshape(nu,nv)
		self.zeta          = zeta.reshape(nu,nv)
		self.phi           = phi.reshape(nu,nv)
		self.rreal         = rreal.reshape(nu,nv)
		self.zreal         = zreal.reshape(nu,nv)
		self.Nx            = Nx.reshape(nu,nv)
		self.Ny            = Ny.reshape(nu,nv)
		self.Nz            = Nz.reshape(nu,nv)
		self.bnreal        = bnreal.reshape(nu,nv)
		self.bcreal        = bcreal.reshape(nu,nv)
		self.bnormal_total = bnormal_total.reshape(nu,nv)
		

	def read_bnorm_harm(self,ext):
		"""Reads a BNORM_HARM file

		This routine reads and initilizes the BNORM class
		with variable information from a BNORM_HARM file.

		Parameters
		----------
		filename : str
			Path to bnorm_harm file.
		"""
		import numpy as np
		self.ext = ext
		f = open("bnorm_harm."+ext,'r')
		lines = f.readlines()
		f.close()
		# First line contains the number of modes
		self.mnmax = int(lines[0])
		lines.pop(0)
		self.m = np.zeros((self.mnmax,1))
		self.n = np.zeros((self.mnmax,1))
		self.bmnc = np.zeros((self.mnmax,1))
		self.bmns = np.zeros((self.mnmax,1))
		for line in lines:
			(txt0,txt1,txt2,txt3,txt4) = line.split()
			mn = int(txt0)-1
			self.m[mn] = int(txt1)
			self.n[mn] = int(txt2)
			self.bmnc[mn] = float(txt3)
			self.bmns[mn] = float(txt4)

	def plot_bnorm_real_histogram(self,ax=None, lsave=False):
		"""Plots a histogram of the Bnormal total values

		This routine plots a histogram of the bnormal total values.

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		ax.hist(self.bnormal_total.flatten(),bins=100)
		ax.set_xlabel(r'$B_{normal}$ [T]')
		ax.set_ylabel('Counts')
		ax.set_title(rf'Bnorm_real Total Values')
		if lsave: pyplot.savefig(self.ext+'/bnorm_real_histogram.png')
		if lplotnow: pyplot.show()

	def print_bnorm_real_stats(self):
		"""Prints some statistics of the Bnormal total values

		This routine prints some statistics of the bnormal total values.
		"""
		import numpy as np
		print(f'BNORMAL REAL STATSITICS:')
		print(f'Minimum: {np.min(self.bnormal_total)}')
		print(f'Maximum: {np.max(self.bnormal_total)}')
		print(f'Mean: {np.mean(self.bnormal_total)}')
		print(f'Standard Deviation: {np.std(self.bnormal_total)}')
		print(f'Variance: {np.var(self.bnormal_total)}')
		# Percentiles
		print(f'25th Percentile: {np.percentile(self.bnormal_total,25)}')
		print(f'50th Percentile: {np.percentile(self.bnormal_total,50)}')
		print(f'75th Percentile: {np.percentile(self.bnormal_total,75)}')

	def plot_bnorm_harm(self,ax=None,type='sin', lsave=False):
		"""Plots the Bnormal spectrum

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		# Plot the spectrum as a tile plot. Each m,n mode is a tile, colored by the value of bmns
		mmax = int(max(np.squeeze(self.m)))
		nmax = int(max(np.squeeze(self.n)))
		bmn = np.zeros((mmax+1,2*nmax+1))
		for mn in range(self.mnmax):
			# m = int(self.m[mn])
			# n = int(self.n[mn]) + nmax
			m = int(self.m[mn].item())
			n = int(self.n[mn].item()) + nmax
			if type == 'sin':
				# bmn[m,n] = self.bmns[mn]
				bmn[m,n] = self.bmns[mn].item()
			elif type == 'cos':
				# bmn[m,n] = self.bmnc[mn]
				bmn[m,n] = self.bmnc[mn].item()
		x = np.linspace(0,mmax,mmax+1)
		y = np.linspace(-nmax,nmax,2*nmax+1)
		# Include colorbar, quadmesh without interpolation
		# The colormap should be symmetric about zero, so that positive and negative values are colored differently
		# Define the norm that ensures the colormap is symmetric about zero
		import matplotlib
		norm = matplotlib.colors.CenteredNorm()
		quadmesh=ax.pcolormesh(x,y,bmn.T,cmap='seismic',norm=norm)
		pyplot.colorbar(quadmesh,label='$B_{normal}$ [arb]',ax=ax)

		ax.set_xlabel('Poloidal Modes (m)')
		ax.set_ylabel('Toroidal Modes (n)')
		ax.set_title(rf'BNORM Harmonic Spectrum (Sin)')
		# if lsave: pyplot.savefig(f'{self.ext}/bnorm_harm_{type}.png')
		if lsave:
			pyplot.savefig(f'{self.ext}/bnorm_harm_{type}.png')
			pyplot.close()
		elif lplotnow:
			pyplot.show()

	def plotBnmnSpectrum(self,ax=None,cmap='jet'):
		"""Plots the Bnormal spectrum for a surface

		This routine plots the bnormal spectrum for a given
		surface.

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		cmap : string (optional)
			Matplotlib colormap.

		Returns
		-------
		quadmesh : matplotlib.collections.Quadmesh
			Quadmesh as produced by pcolormesh
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		# Array extents
		mmax = int(max(np.squeeze(self.xm)))
		nmax = int(max(np.squeeze(self.xn)))
		# Sort BMN into array
		bmn = np.zeros((mmax+1,2*nmax+1))
		for mn in range(self.mnmax):
			m = int(self.xm[mn])
			n = int(self.xn[mn]) + nmax
			bmn[m,n] = self.bnmns[0,mn]
		#Plot
		x = np.linspace(0,mmax,mmax+1)
		y = np.linspace(-nmax,nmax,2*nmax+1)
		quadmesh=ax.pcolormesh(x,y,np.log10(np.abs(bmn.T)),cmap=cmap,shading='gouraud')
		ax.set_xlabel('Poloidal Modes (m)')
		ax.set_ylabel('Toroidal Modes (n)')
		ax.set_title(rf'BNORM Normal Field')
		pyplot.colorbar(quadmesh,label='$log_{10}$[arb]',ax=ax)
		if lplotnow: pyplot.show()
		return quadmesh

	def plot_bnorm_real_total(self,ax=None):
		"""Plots the Bnormal spectrum for a surface

		This routine plots the bnormal spectrum for a given
		surface.

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		x = self.theta[:,0]
		y = self.phi[0,:]
		hmesh=ax.pcolormesh(y,x,self.bnormal_total.T)
		ax.set_xlabel('Toroidal Angle (phi) [rad]')
		ax.set_ylabel('Poloidal Angle (phi) [rad]')
		ax.set_title(rf'Total B-Normal Field')
		# pyplot.colorbar(hmesh,label='$log_{10}$[arb]',ax=ax)
		if lplotnow: pyplot.show()

	def plotBsurf(self,ax=None,cmap='jet',lsave=False):
		"""Plots the Bnormal on a surface

		This routine plots the bnormal.

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		cmap : string (optional)
			Matplotlib colormap.

		Returns
		-------
		quadmesh : matplotlib.collections.Quadmesh
			Quadmesh as produced by pcolormesh
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True
		theta = np.deg2rad(np.linspace([0],[360],360))
		zeta  = np.deg2rad(np.linspace([0],[360],256))
		b = self.sfunct(theta,zeta,self.bnmns,self.xm,self.xn)
		quadmesh=ax.pcolormesh(np.squeeze(zeta),np.squeeze(theta),np.squeeze(b[0,:,:]),cmap=cmap,shading='gouraud')
		ax.set_xlabel(r'Toroidal Angle ($\zeta=\phi N_{fp}$) [rad]')
		ax.set_ylabel(r'Poloidal Angle ($\theta$) [rad]')
		ax.set_title(rf'BNORM')
		pyplot.colorbar(quadmesh,label=r'$B_{normal}$ [arb]',ax=ax)
		if lsave: pyplot.savefig(f'{self.ext}/bnorm_surf.png')
		if lplotnow: pyplot.show()
		return quadmesh

	def plotBrealsurf(self,ax=None,cmap='jet',lsave=False):
		"""Plots the Bnormal on a surface

		This routine plots the bnormal.

		Parameters
		----------
		ax : axes (optional)
			Matplotlib axes object to plot to.
		cmap : string (optional)
			Matplotlib colormap.

		Returns
		-------
		scatter : matplotlib.collections.PathCollection
			PathCollection as produced by scatter
		"""
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			ax = pyplot.axes()
			lplotnow = True

		scatter = ax.scatter(self.theta.flatten(), self.phi.flatten(), c=self.bnormal_total.flatten(), cmap=cmap, alpha=0.7)
		ax.set_xlabel(r'Poloidal Angle ($\theta$) [rad]')
		ax.set_ylabel(r'Toroidal Angle ($\phi$) [rad]')
		ax.set_title(rf'BNORM Real')
		pyplot.colorbar(scatter,label=r'$B_{normal}$ [T]',ax=ax)
		
		if lsave: pyplot.savefig(f'{self.ext}/bnorm_real_surf.png')
		if lplotnow: pyplot.show()
		return scatter

	def plotBrealsurf_3D(self,ax=None,cmap='jet',lsave=False):
		# Scatters the bnormal values on the surface, colored by the bnormal value
		import numpy as np
		import matplotlib.pyplot as pyplot
		lplotnow = False
		if not ax:
			# Create 3D axes
			ax = pyplot.axes(projection='3d')
			lplotnow = True

		print(self.xreal.shape, self.yreal.shape, self.zzreal.shape, self.full_bnormal_total.shape)
		ax.scatter(self.xreal, self.yreal, self.zzreal, c=self.full_bnormal_total, cmap='viridis', alpha=0.7)
		# Colorbar
		pyplot.colorbar(ax.collections[0], label='B_n (T)', ax=ax)
		ax.set_aspect('equal')

		if lsave: pyplot.savefig(f'{self.ext}/bnorm_real_surf.png')
		if lplotnow: pyplot.show()
		return ax

# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)









