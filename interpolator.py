from scipy import signal


from pylab import *


fig,(ax1,ax2)=subplots(2,1,num=1,clear=True)

n=420
u=11
d=6


if 1:
    s=sin(2*pi*arange(n)/30)+sin(2*pi*arange(n)/25)
else:
    s=zeros(420)
    s[10]=1

ax2.plot(arange(len(s)),s,'.-')

#______ calculate equivalent filter
far=lambda mu:array([mu*(-0.5+mu/2),-mu/2-mu**2/2+1,mu*(1.5-mu/2),-mu/2+mu**2/2])
far=lambda mu:array([1-mu,mu])

#______ direct approach
y0=[]
lenfar=len(far(0))
for i in range(len(s)*u//d-10):
    j=i/u*d
    y0.append(s[int(j):int(j)+lenfar]@far(j%1))
    
ax2.plot(arange(len(y0))*d/u+1,real(y0),'.-');


#_____ upfirdn approach
fir=reshape(list(zip(*[far(i/u) for i in range(u,0,-1)])),[-1])
y1=signal.upfirdn(fir,s,u,d)


#_____ Frequency response
FIR=fft(fftshift(concatenate([zeros(n*u//2-len(fir)//2),fir,zeros(n*u//2-len(fir)//2)])))
FIRds=real(concatenate([FIR[-n*u//d//2:],FIR[:n*u//d//2]]))/u
f=linspace(-u/d/2,u/d/2,n*u//d)
ax1.plot(f,FIRds)


#_____ FFT approach
S=tile(fft(s),u)
S*=FIR
SS=sum(reshape(S,[d,-1]),0)
FIRacc=sum(reshape(FIR,[d,-1]),0)


ax1.plot(linspace(-u/d/2,u/d/2,len(FIRacc)),real(fftshift(FIRacc))/u)
#ax1.plot(linspace(-u/2,u/2,len(FIR)),real(fftshift(FIR)))
ax1.grid()
ax1.legend(['upsample only','with ds overlap'])

y2=ifft(SS)/d
#plot(20*log10(abs(fft(fir,1024))+.01))


ax2.plot(arange(len(y1))*d/u-u+d,y1,'.-');
ax2.plot(arange(len(y2))*d/u-u+d,real(y2),'.-');
ax2.legend(['original','direct','upFIRn','FFT'])