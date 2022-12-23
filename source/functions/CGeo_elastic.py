import numpy as np
import matplotlib.pyplot as plt
# The Hilbert transform is used to create the wavelet
from scipy.signal import hilbert

"""
	This program simulates wave propagation in elastic media
	It has four classes with their own data and functions:
	
	Source: simulates the wave source
	Derivatives: Solves derivatives by finite differences
	Elastic_model: Sets the elastic model and makes plots
	Elastic_waves: Simulates the wave propagation
	
	Theory and numerical solution from Virieux (1986)
	
	Written by Wiktor Weibull (wiktor.w.weibull@uis.no)
"""

class Source:
	def __init__(self, nt, dt, sx, sz):
		self.nt = nt
		self.dt = dt
		self.sx = sx
		self.sz = sz
		self.t = np.linspace(0, (nt-1)*dt, nt)
		self.wav = np.zeros([nt])
	
	def Ricker(self, f0, t0, phase):
		arg = np.pi*np.pi*f0*f0*(self.t-t0)**2
		self.wav[:] = (1-2*arg)*np.exp(-arg)
		# Rotate phase
		self.wav = np.real(np.fft.ifft(np.fft.fft(hilbert(self.wav)*(np.exp(1j*phase)))))
		return self.wav, self.t
	
	def plot(self):
		fig, ax = plt.subplots(figsize=(8,6))
		ax.plot(self.t,self.wav)
		ax.set_xlabel("Time (s)")
		ax.set_ylabel("Amplitude") 
		ax.set_title("Wavelet")
		ax.grid()
		return fig, ax
	
class Derivatives:
	def __init__(self): 
		pass
	
	# Compute discrete forward derivatives
	def Dxfw(self,M, dx):
		nz,nx = M.shape
		D = np.zeros(M.shape)
		D[0:nz,0:nx-1] = M[0:nz,1:nx] - M[0:nz,0:nx-1]
		D = D/(dx)
		return D
	
	def Dzfw(self,M, dz):
		nz,nx = M.shape
		D = np.zeros(M.shape)
		D[0:nz-1,0:nx] = M[1:nz,0:nx] - M[0:nz-1,0:nx]
		D = D/(dz)
		return D
	
	# Compute discrete backward derivatives
	def Dxbw(self,M, dx):
		nz,nx = M.shape
		D = np.zeros(M.shape)
		D[0:nz,1:nx] = M[0:nz,1:nx] - M[0:nz,0:nx-1]
		D = D/(dx)
		return D
	
	def Dzbw(self,M, dz):
		nz,nx = M.shape
		D = np.zeros(M.shape)
		D[1:nz,0:nx] = M[1:nz,0:nx] - M[0:nz-1,0:nx]
		D = D/(dz)
		return D
	
class Elastic_model:
	def __init__(self, vp, vs, rho, dx, dz, ox, oz):
		self.vp = vp
		self.vs = vs
		self.rho = rho
		nz, nx = vp.shape
		self.nz = nz
		self.nx = nx
		self.dx = dx
		self.dz = dz
		self.ox = ox
		self.oz = oz
		
		for i in range(0,nx):
			for j in range (0,nz):
				if(vp[j,i] < vs[j,i]*np.sqrt(4./3.)):
					raise Exception("This elastic model is not physically possible" 
					"since the Possion ratio is less than 0.5.")
		
		# Computing moduli
		self.L2M = rho*vp*vp
		self.M = rho*vs*vs
		self.L = self.L2M - 2*self.M;

		# Staggering
		self.B_x = 1./rho
		self.B_x[:,0:nx-1] = 1./(0.5*(rho[:,0:nx-1]+rho[:,1:nx]))
		self.B_z = 1./rho
		self.B_z[0:nz-1,:] = 1./(0.5*(rho[0:nz-1,:]+rho[1:nz,:]))
		self.Mxz = self.M
		self.Mxz[:,0:nx-1] = 0.5*(self.M[:,0:nx-1] + self.M[:,1:nx])
		self.Mxz[0:nz-1,:] = 0.5*(self.Mxz[0:nz-1,:] + self.Mxz[1:nz,:])
		
		# Adding free surface condition at the top
		self.B_x[1,:] = 2.0*self.B_x[1,:]
		self.L[1,:] = 0.0
		self.L2M[1,:] = self.M[1,:]
	
	def plot(self):
		x0 = self.ox
		xend = x0 + (self.nx-1)*self.dx
		z0 = self.oz
		zend = z0 + (self.nz-1)*self.dz
		fig, ax = plt.subplots(1,3,sharey=True,figsize=(15,7.5))
		im0 = ax[0].imshow(self.vp, extent=[x0, xend, zend, z0], aspect=2, cmap="rainbow")
		ax[0].set_xlabel("X (m)")
		ax[0].set_ylabel("Z (m)")
		ax[0].set_title("P-wave velocity");
		cbar0=fig.colorbar(im0, ax=ax[0], orientation="horizontal")
		cbar0.ax.set_xlabel("Vp (m/s)")
		im1 = ax[1].imshow(self.vs, extent=[x0, xend, zend, z0], aspect=2, cmap="rainbow")
		ax[1].set_xlabel("X (m)")
		ax[1].set_title("S-wave velocity");
		cbar1=fig.colorbar(im1, ax=ax[1], orientation="horizontal")
		cbar1.ax.set_xlabel("Vs (m/s)")
		im2 = ax[2].imshow(self.rho, extent=[x0, xend, zend, z0], aspect=2, cmap="rainbow")
		ax[2].set_xlabel("X (m)")
		ax[2].set_title("Density");
		cbar2=fig.colorbar(im2, ax=ax[2], orientation="horizontal")
		cbar2.ax.set_xlabel(r"rho (kg/m$^3$)")
		return fig, ax

class Elastic_waves:
	def __init__(self, Model, nt, dt):
		self.nz = Model.nz
		self.nx = Model.nx
		self.dx = Model.dx
		self.dz = Model.dz
		self.ox = Model.ox
		self.oz = Model.oz
		
		self.Sxx = np.zeros([self.nz,self.nx]) # Stress component
		self.Szz = np.zeros([self.nz,self.nx]) # Stress component
		self.Sxz = np.zeros([self.nz,self.nx]) # Stress component
		self.Vx = np.zeros([self.nz,self.nx]) # Particle velocity component
		self.Vz = np.zeros([self.nz,self.nx]) # Particle velocity component
		self.nt=nt
		self.dt=dt
	
	def Courant_stability(self, vpmax):
		dtstab = np.sqrt(self.dx**2 + self.dz**2)/(np.sqrt(2)*vpmax)
		return dtstab
	
	def insertPressure(self, source, it):
		sx = source.sx
		sz = source.sz
		wav = source.wav[it]
		self.Sxx[sz,sx] = self.Sxx[sz,sx] + self.dt*wav
		self.Szz[sz,sx] = self.Szz[sz,sx] + self.dt*wav
	
	def insertForce(self, source, model, direction, it):
		sx = source.sx
		sz = source.sz
		wav = source.wav[it]
		dx = np.sin(direction*np.pi/180.)
		dz = -np.cos(direction*np.pi/180.)
		Bx = model.B_x
		Bz = model.B_z
		self.Vx[sz,sx-1] = self.Vx[sz,sx-1] + dx*0.5*self.dt*Bx[sz,sx-1]*wav
		self.Vx[sz,sx] = self.Vx[sz,sx] + dx*0.5*self.dt*Bx[sz,sx]*wav
		self.Vz[sz-1,sx] = self.Vz[sz-1,sx] + dz*0.5*self.dt*Bz[sz-1,sx]*wav
		self.Vz[sz,sx] = self.Vz[sz,sx] + dz*0.5*self.dt*Bz[sz,sx]*wav
	
	def insertStress(self, source, which, it):
		sx = source.sx
		sz = source.sz
		wav = source.wav[it]
		if(which == "xx"):
			self.Sxx[sz,sx] = self.Sxx[sz,sx] + self.dt*wav
		if(which == "zz"):
			self.Szz[sz,sx] = self.Szz[sz,sx] + self.dt*wav
		if(which == "xz"):
			self.Sxz[sz,sx] = self.Sxz[sz,sx] + 0.25*self.dt*wav
			self.Sxz[sz-1,sx] = self.Sxz[sz-1,sx] + 0.25*self.dt*wav
			self.Sxz[sz,sx-1] = self.Sxz[sz,sx-1] + 0.25*self.dt*wav
			self.Sxz[sz-1,sx-1] = self.Sxz[sz-1,sx-1] + 0.25*self.dt*wav
	
	def recordPressure(self, rz):
		return (self.Sxx[rz,:] + self.Szz[rz,:])
	
	def recordVelocity(self, rz, direction):
		record = np.zeros([self.nx])
		if(direction == "x"):
			record[1:self.nx] = 0.5*self.Vx[rz,0:self.nx-1]+0.5*self.Vx[rz,1:self.nx]
			return(record)
		if(direction == "z"):
			record[0:self.nx] = 0.5*self.Vz[rz-1,:]+0.5*self.Vz[rz,:]
			return(record)
	
	def forwardStep(self,Derivative,Model):
		self.forwardstepVelocity(Derivative,Model)
		self.forwardstepStress(Derivative,Model)
	
	# Solve elastodynamic equations
	def forwardstepVelocity(self,Derivative,Model):
		dx = self.dx
		dz = self.dz
		dt = self.dt
		nx = self.nx
		nz = self.nz
		
		## Forward step Vx
		der = Derivative.Dxfw(self.Sxx,dx)
		self.Vx = self.Vx + dt*Model.B_x*der
		
		der = Derivative.Dzbw(self.Sxz,dz) 
		self.Vx = self.Vx + dt*Model.B_x*der
		
		## Forward step Vz
		der = Derivative.Dxbw(self.Sxz,dx)
		self.Vz = self.Vz + dt*Model.B_z*der
		
		der = Derivative.Dzfw(self.Szz,dz)
		self.Vz = self.Vz + dt*Model.B_z*der
	
	# Solve elastodynamic equations
	def forwardstepStress(self,Derivative,Model):
		dx = self.dx
		dz = self.dz
		dt = self.dt
		nx = self.nx
		nz = self.nz
		
		# Forward step normal stresses
		der = Derivative.Dxbw(self.Vx,dx)
		self.Sxx = self.Sxx + dt*Model.L2M*der
		self.Szz = self.Szz + dt*Model.L*der
		
		der = Derivative.Dzbw(self.Vz,dz)
		self.Sxx = self.Sxx + dt*Model.L*der
		self.Szz = self.Szz + dt*Model.L2M*der
		
		# Free surface condition at the top
		self.Szz[1,:] = 0.0
		
		# Forward step shear stress
		der = Derivative.Dxfw(self.Vz,dx)
		self.Sxz = self.Sxz + dt*Model.Mxz*der
		
		der = Derivative.Dzfw(self.Vx,dz)
		self.Sxz = self.Sxz + dt*Model.Mxz*der
