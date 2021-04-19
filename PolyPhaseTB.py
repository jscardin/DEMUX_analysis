#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 10:25:27 2021

@author: cardinal
"""
from pylab import *
from quicktions import Fraction


fig,[ax1,ax2]=subplots(2,1,num=2,clear=True)

far=lambda mu:array([mu*(-0.5+mu/2),-mu/2-mu**2/2+1,mu*(1.5-mu/2),-mu/2+mu**2/2]) #Farrow interpolation
#far=lambda mu:array([1-mu,mu])   # Linear interpolation

def interp_farrow(xx,yy):
    y=[]
    yy=concatenate([[yy[0]],yy,[yy[-1]]*2])  # to extrapolate
    for x in xx:
        y.append(far(x%1)@yy[int(x):int(x)+4])
    return(array(y))


PolyWidth=4
FIRlen=32*PolyWidth
OS=2.6
OSfin=4
n=1000

(up,dn)=Fraction(OSfin/OS).limit_denominator().as_integer_ratio()
print(up,dn)

a=real(fftshift(ifft(RRC1(FIRlen,FIRlen,1/OS/PolyWidth,0.2))))
#b=interp(arange(len(a)*up//PolyWidth)/up*PolyWidth,arange(len(a)),a)
b=interp_farrow(arange(len(a)*up//PolyWidth)/up*PolyWidth,a)
c=concatenate([zeros(n//2-len(b)//2),b,zeros(n//2-len(b)//2)])


#plot(b[::dn],'.-')
ax1.plot(arange(len(a))/PolyWidth,a,'.-',lw=1)
ax1.plot(arange(len(b[::dn]))/len(b)*dn*FIRlen/PolyWidth,b[::dn],'.-')
ax2.plot(abs(fft(b)))
ax2.plot(abs(fft(fftshift(c))))
