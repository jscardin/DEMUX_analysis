
from pylab import *
from scipy.signal import upfirdn

N=20
u=4
d=2
FIR=array([0.25,0.5,0.75])
FIR=concatenate([FIR,[1],flipud(FIR)])

t=arange(N)
tt=linspace(0,N,N*u//d,endpoint=False)
x=sin(2*pi*t/10)
y=upfirdn(FIR,x,u,d,mode="wrap")

fir=figure(1,clear=True)
plot(t,x,'.-')
plot(tt,y[:len(tt)],'.-')

