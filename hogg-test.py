#!/usr/bin/env python
#Create Hogg plots for demonstration purposes

from cosmos import *
from pylab import *

z=arange(0.1,5,0.1)
a=Cosmo(omega0=1.0,omegalambda=0.0,verbose=0)
b=Cosmo(omega0=0.05,omegalambda=0.0,verbose=0)
c=Cosmo(omega0=0.2,omegalambda=0.8,verbose=0)

dh1,th1=a.dhubble(),a.thubble()
dh2,th2=b.dhubble(),b.thubble()
dh3,th3=c.dhubble(),c.thubble()

#plot 1: z vs dcomovingtransverse
ioff()
y1,y2,y3=[],[],[]
for zz in z:
 y1.append(a.dcomovingtransverse(zz)/dh1)
 y2.append(b.dcomovingtransverse(zz)/dh2)
 y3.append(c.dcomovingtransverse(zz)/dh3)
ax=subplot(111)
ax.plot(z,y1,'--',z,y2,'--',z,y3,'--')
xlim(0,5)
ylim(0,3)

#majorLocator=MultipleLocator(1)
#majorFormatter=FormatStrFormatter('%.2f')
#minorLocator=MultipleLocator(.1)
#ax.xaxis.set_major_locator(majorLocator)
#ax.xaxis.set_major_formatter(majorFormatter)
#ax.xaxis.set_minor_locator(minorLocator)
#ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
#ax.yaxis.set_minor_locator(minorLocator)

title('Plot 1')
xlabel('redshift z')
ylabel('Proper Motion Distance D_dM/D_dH')
#legend(('1','2','3'),loc='upper left')
savefig('pyhogg_plot1.png')
savefig('pyhogg_plot1.ps')

#plot 2: z vs dangular/dh
clf()
y1,y2,y3=[],[],[]
for zz in z:
 y1.append(a.dangular(zz)/dh1)
 y2.append(b.dangular(zz)/dh2)
 y3.append(c.dangular(zz)/dh3)
ax=subplot(111)
ax.plot(z,y1,'--',z,y2,'--',z,y3,'--')
xlim(0,5)
#ylim(0,3)
title('Plot 2')
xlabel('redshift z')
ylabel('Angular Diameter Distance D_dA_n/D_dH')
#legend(('1','2','3'),loc='upper left')
savefig('pyhogg_plot2.png')
savefig('pyhogg_plot2.ps')

#plot 3: z vs dluminosity/dh
clf()
y1,y2,y3=[],[],[]
for zz in z:
 y1.append(a.dluminosity(zz)/dh1)
 y2.append(b.dluminosity(zz)/dh2)
 y3.append(c.dluminosity(zz)/dh3)
ax=subplot(111)
ax.plot(z,y1,'--',z,y2,'--',z,y3,'--')
xlim(0,5)
#ylim(0,3)
title('Plot 3')
xlabel('redshift z')
ylabel('Luminosity Distance D_dL_n/D_dH')
#legend(('1','2','3'),loc='upper left')
savefig('pyhogg_plot3.png')
savefig('pyhogg_plot3.ps')

#plot 4: z vs dmodulus(z) +5*log10(red100) ??
clf()
y1,y2,y3=[],[],[]
for zz in z:
 y1.append(a.dmodulus(zz)+5*log10(a.h100))
 y2.append(b.dmodulus(zz)+5*log10(a.h100))
 y3.append(c.dmodulus(zz)+5*log10(a.h100))
ax=subplot(111)
ax.plot(z,y1,'--',z,y2,'--',z,y3,'--')
xlim(0,5)
#ylim(0,3)
title('Plot 4')
xlabel('redshift z')
ylabel('Comoving Volume Element ... [stuff]')
#legend(('1','2','3'),loc='upper left')
savefig('pyhogg_plot4.png')
savefig('pyhogg_plot4.ps')

#plot 5: z vs dvcomoving(z)/(dhubble()**3)
clf()
y1,y2,y3=[],[],[]
for zz in z:
 y1.append(a.dvcomoving(zz)/dh1**3)
 y2.append(b.dvcomoving(zz)/dh2**3)
 y3.append(c.dvcomoving(zz)/dh3**3)
ax=subplot(111)
ax.plot(z,y1,'--',z,y2,'--',z,y3,'--')
xlim(0,5)
#ylim(0,3)
title('Plot 5')
xlabel('redshift z')
ylabel('Comoving Volume Element ... [stuff]')
#legend(('1','2','3'),loc='upper left')
savefig('pyhogg_plot5.png')
savefig('pyhogg_plot5.ps')

#plot 6: ??
clf()
y1,y2,y3,y4,y5,y6=[],[],[],[],[],[]
for zz in z:
 cage1=a.getage(zz)/th1
 cage2=b.getage(zz)/th2
 cage3=c.getage(zz)/th3
 clook1=a.getage(0)/th1-cage1
 clook2=b.getage(0)/th2-cage2
 clook3=c.getage(0)/th2-cage3
 y1.append(cage1)
 y2.append(cage2)
 y3.append(cage3)
 y4.append(clook1)
 y5.append(clook2)
 y6.append(clook3)
ax=subplot(111)
ax.plot(z,y1,',',z,y2,'.',z,y3,'--',z,y4,',',z,y5,'.',z,y6,'--')
xlim(0,5)
#ylim(0,3)
title('Plot 6')
xlabel('redshift z')
ylabel('Lookback Time ... [stuff]')
#legend(('1','2','3'),loc='upper left')
savefig('pyhogg_plot6.png')
savefig('pyhogg_plot6.ps')

#plot 7: z vs (1+z)**2/epeebles(z)
clf()
y1,y2,y3=[],[],[]
for zz in z:
 y1.append((1+zz)**2/a.epeebles(zz))
 y2.append((1+zz)**2/b.epeebles(zz))
 y3.append((1+zz)**2/c.epeebles(zz))
ax=subplot(111)
ax.plot(z,y1,'--',z,y2,'--',z,y3,'--')
xlim(0,5)
#ylim(0,3)
title('Plot 7')
xlabel('redshift z')
ylabel('Dimensionless Intersection Probability dP/dz')
#legend(('1','2','3'),loc='upper left')
savefig('pyhogg_plot7.png')
savefig('pyhogg_plot7.ps')
