"""Kick-ass Python Cosmology package meant to mimic the IDL cosmology
    package of Leonidas Moustakas and John Moustakas.
"""

import numpy
from math import *
from pygsl import integrate as integ
from pygsl import interpolation as interp
from pygsl import spline

#general functions
def convert_distance(dunit='pc'):
 """ Returns distance conversion factor.  Default unit is pc, corresponding to
     a scale factor of 1.  Other units are cm, m, kpc, Mpc."""
 if dunit is 'pc':    convert=1.0  #default unit is pc
 elif dunit is 'cm':  convert=3.086e18
 elif dunit is 'm':   convert=3.086e16
 elif dunit is 'kpc': convert=1.0e-3
 elif dunit is 'Mpc': convert=1.0e-6
 return convert

def convert_time(tunit='Gyr'):
 """ Returns time conversion factor.  Default unit is Gyr, corresponding to
     a scale factor of 1.  Other units are s, yr, Myr."""
 secyr=365.25*24*3600
 if tunit is 'Gyr':   convert=1.0 #default unit is Gyr
 elif tunit is 's':   convert=1.0e9*secyr
 elif tunit is 'yr':  convert=1.0e9
 elif tunit is 'Myr': convert=1.0e3
 return convert

def asinh(x):
 """ The inverse hyperbolic sine function."""
 return log(x + sqrt(1.+x*x))

def cube(x):
 """ Returns the cube of x."""
 x=float(x)
 return x**3

def cuberoot(x):
 """Returns the cube root of x"""
 return x**(1/3.)

def sqr(x):
 """ Returns the square of x."""
 x=float(x)
 return x**2

def romberg(f, a, b, eps = 1E-8):
    """Approximate the definite integral of f from a to b by Romberg's method.
    eps is the desired accuracy.

    romberg is the equivalent of the IDL function qromb."""
    R = [[0.5 * (b - a) * (f(a) + f(b))]]  # R[0][0]
    n = 1
    while True:
        h = float(b-a)/2**n
        R.append((n+1)*[None])  # Add an empty row.
        R[n][0] = 0.5*R[n-1][0] + h*sum(f(a+(2*k-1)*h) for k in range(1, 2**(n-1)+1)) # for proper limits
        for m in range(1, n+1):
            R[n][m] = R[n][m-1] + (R[n][m-1] - R[n-1][m-1]) / (4**m - 1)
        if abs(R[n][n-1] - R[n][n]) < eps:
            return R[n][n]
        n += 1

class Cosmo:
 """Kick-ass Python Cosmology Object developed to mimic the IDL cosmology
    package of Leonidas Moustakas and John Moustakas.  All functions appear 
    in David Hogg's paper at http://arxiv.org/abs/astro-ph/9905116v4.  The
    accompanying hogg-test.py script reproduces plots from this paper.

 Cosmo is initialized with the following parameters:
   omega0=Density of matter; default is 0.3
   omegalambda=Density of cosmological constant; default is 0.7
   h100=0.01*Hubble constant; default is 0.7
   dunit = Distance unit; default is 'pc'
   tunit = Time unit; default is 'Gyr'
   verbose = Verbose mode; default is True.  Turn off with False.

 Make the Default cosmology:
   >>>a=Cosmo()
   Omega Matter    = 0.300000
   Omega Lambda    = 0.700000
   H_0 / 100       = 0.700000

 Make a wacky custom cosmology:
   >>>b=Cosmo(omega0=0.654,omegalambda=0.456,h100=0.42)
   Omega Matter    = 0.654000
   Omega Lambda    = 0.456000
   H_0 / 100       = 0.420000


 Changing units:
  Distance and time units are handled with dunit and tunit.  Units can be
  specificied both in the Cosmo constructor and dyamically at any later time.
  Setting dunit to any one of 'cm', 'm', 'pc', 'kpc', or 'Mpc' specifies the 
  output unit of all subsequent distance calculations until dunit is changed 
  again.  Setting tunit to any one of 's', 'yr', 'Myr', or 'Gyr' will do 
  likewise for all time calculations.

  Once a Cosmo() object has been specified, the units can be changed by 
  adjusting the dunit of that object directly:
   >>>a=Cosmo(verbose=0)
   >>>a.dhubble()
   4282749399.9999995
   >>>a.dunit='Mpc'
   >>>a.dhubble()
   4282.7493999999997
 """

 def __init__(self,omega0=0.3,omegalambda=0.7,h100=0.7,verbose=1,dunit='pc',tunit='Gyr'):
  """Cosmo object constructor.  Define a Cosmo object with the default or a
     custom cosmology, e.g.:
      >>>>a=Cosmo()
      >>>>b=Cosmos(omega0=0.654,omegalambda=0.456,h100=0.42)
  """
  self.omega0=omega0
  self.omegalambda=omegalambda
  self.h100=h100
  self.verbose=verbose
  self.dunit=dunit
  self.tunit=tunit
  if verbose:
   print 'Omega Matter\t= %f\nOmega Lambda\t= %f\nH_0 / 100\t= %f'%(omega0,omegalambda,h100)

 def update(self,omega0=0.3000,omegalambda=0.7,h100=0.70,verbose=1,dunit='pc',tunit='Gyr'):
  """Update multiple Cosmo properties with a single function call.
  """
  self.omega0=omega0
  self.omegalambda=omegalambda
  self.h100=h100
  self.verbose=verbose
  self.dunit=dunit
  self.tunit=tunit
  if verbose:
   print 'Omega Matter\t= %f\nOmega Lambda\t= %f\nH_0 / 100\t= %f'%(omega0,omegalambda,h100)
   print 'Distance unit\t= %s\nTemporal unit\t= %s'%(dunit,tunit)

 def spill(self):
  """Print the properties of a Cosmo object all at once.  Also returns a list
     of the properties: [omega0,omegalambda,h100,dunit,tunit]
  """
  print 'Omega Matter\t= %f\nOmega Lambda\t= %f\nH_0 / 100\t= %f'%(self.omega0,self.omegalambda,self.h100)
  print 'Distance unit\t= %s\nTemporal unit\t= %s'%(dunit,tunit)
  return [self.omega0,self.omegalambda,self.h100,self.dunit,self.tunit]

 def epeebles(self,z): 
  """A convenience function for internal use having little external value.
     Equation 14 in the Hogg paper.  
  """
  omegaR = 1. - self.omega0-self.omegalambda
  ez = sqrt(self.omega0*cube(1.+z)+omegaR*sqr(1.+z)+self.omegalambda)
  return ez

 def dhubble(self): 
  """The Hubble distance, defined as c/H_0.  The units of the returned values
     correspond to dunit.
  """
  c=2.99792458e5 #(km/s)
  dh = c/self.h100/100. # (Mpc)
  conv=convert_distance(self.dunit) # (pc)
  return conv*dh*1.0e6 # [1]*[Mpc]*[10e6 pc/1 Mpc] = pc

 def thubble(self): 
  """The Hubble time, defined as 1/H_0.  The units of the returned values
     correspond to tunit.
  """
  #recall, 1./(1 km/s/mpc) = 9.77813 Gyr
  conv=convert_time(self.tunit) #(Gyr)
  return conv*9.77813/self.h100

 def dangular(self,z): #returns dunit per radian at z
  """Angular diameter distance, defined as ratio of physical transverse size
     to angular size in radians as a function of z.  Returns dunit per radian
     at redshift z.  Useful for calculating size per arcsecond.  Equation 18 
     in the Hogg paper.
     ...
     1 radian = 206,264.806 arcseconds
  """
  return self.dcomovingtransverse(z)/(1.+z)

 def dangulardiff(self,z1,z2):
  """Angular diameter distance between two objects at redshifts z1 and z2.
     Hogg paper Equation 19.
  """
  omegaR = 1. -self.omega0-self.omegalambda
  if omegaR >= 0:
   dm1=self.dcomovingtransverse(z1)
   dm2=self.dcomovingtransverse(z2)
   dh = self.dhubble()
   da12 = (1./(1.+z2))*( dm2*sqrt(1.+omegaR*sqr(dm1)/sqr(dh)) - dm1*sqrt(1.+omegaR*sqr(dm2)/sqr(dh)) )
   return da12
  else:
   print 'No formula available for a negative curvature density parameter, jabroni'
   return -1.

 def dlosfunc(self,z):
  """Function for internal use.  Returns 1 over epeebles(z).
  """
  return 1./self.epeebles(z)

 def dcomovinglos(self,z):
  """Line-of-sight comoving distance.  Equation 15 in the Hogg paper.
  """
  #romberg is equivalent of idl qromb
  z=float(z)
  return self.dhubble()*romberg(self.dlosfunc,0,z) 

 def dcomovingtransverse(self,z):
  """Transverse co-moving distance.  Equation 16 in the Hogg paper.
  """
  omegaR=1.-self.omega0-self.omegalambda
  if omegaR > 0:
   dct=self.dhubble()/sqrt(omegaR) * sinh(sqrt(omegaR)*self.dcomovinglos(z)/self.dhubble())
  elif omegaR == 0:
   dct=self.dcomovinglos(z)
  elif omegaR < 0:
   dct =self.dhubble()/sqrt(abs(omegaR))*sin( sqrt(abs(omegaR))*self.dcomovinglos(z)/self.dhubble() )
  return dct

 def dluminosity(self,z):
  """Luminosity distance as a function of z.  Equation 20 in the Hogg paper.
  """
  return self.dcomovingtransverse(z)*(1.+z)

 def dmodulus(self,z):
  """Distance modulus as a function of z.  Hogg paper Equation 25.
     Abs Mag = app mag - distance modulus - K correction
  """
  conv=convert_distance(self.dunit)
  #change dunit temporarily
  dunit=self.dunit
  self.dunit='pc'
  result=conv*5.*log10( self.dluminosity(z)/10. )
  self.dunit=dunit
  return result

 def dvcomoving(self,z):
  """Differential comoving volume element as a function of z.  Equation 28
     in the Hogg paper.
  """
  return self.dhubble() * sqr(1.+z) * sqr(self.dangular(z)) / self.epeebles(z)

 def vcomoving(self,z):
  """Integrated comoving volume as a function of z.  Hogg paper Equation 29.
  """
  omegaR=1. - self.omega0 - self.omegalambda
  dh=self.dhubble()
  dm=self.dcomovingtransverse(z)
  if omegaR > 0:
   vc = (4.*pi*cube(dh))/(2.*omegaR) * ( dm/dh*sqrt(1.+omegaR*sqr(dm)/sqr(dh)) - 1./sqrt(abs(omegaR))*asinh(sqrt(abs(omegaR))*dm/dh) )
  elif omegaR == 0:
   vc = (4.*pi/3.)*cube(dm)
  elif omegaR < 0.:
   vc = (4.*pi*cube(dh))/(2.*omegaR) * ( dm/dh*sqrt(1.+omegaR*sqr(dm)/sqr(dh)) - 1./sqrt(abs(omegaR))*asin(sqrt(abs(omegaR))*dm/dh) )
  return vc

 def getage(self,z): 
  """Look-back time as a function of redhshift.  Returns time in tunit.  
     Hogg paper Equation 30.
  """
  def agefunc(z2,y):
   return 1./ (self.epeebles(z2)*(1.+z2))
  
  #gsl baggage
  gsl_agefunc=integ.gsl_function(agefunc,None) #extra pygsl baggage
  w=integ.workspace(1000000)

  try: #check if z is list or array
   output=[]
   for zz in z:
    flag,result,error=integ.qagiu(gsl_agefunc,zz,1e-8,1e-8,100000,w) 
    output.append(self.thubble()*result)
   output=numpy.array(output)
   return output
  except TypeError: #z is single value
   flag,result,error=integ.qagiu(gsl_agefunc,z,1e-8,1e-8,100000,w) 
   output=self.thubble()*result
  return output

 def getredshift(self,age0): #given an age, return redshift via interpolation
  """Returns the look-back time for a given z.  Age is determined by
     interpolation with the GNU Science Library.  Output is consistent with
     getage.
  """
  #outperforms idl version
  z1,z2,dz=0.1,10.0,0.1
  z=numpy.arange(z1,z2,dz)
  z=z.tolist()
  z.reverse() #reverse redshift order so ages are in increasing order
  ages=self.getage(z)
  interpol=interp.cspline( len(ages) )  #define interpolation object
  interpol.init( ages,z )              #initialize with age the dependent variable
  return interpol.eval(age0)
