from scipy import signal


from pylab import *


fig,(ax1,ax2)=subplots(2,1,num=1,clear=True)

u=3
d=2
far=lambda mu:array([mu*(-0.5+mu/2),-mu/2-mu**2/2+1,mu*(1.5-mu/2),-mu/2+mu**2/2])

fir=reshape(list(zip(*[far(i/u) for i in range(u,0,-1)])),[-1])

if 1:
    s=sin(2*pi*arange(420)/30)+sin(2*pi*arange(420)/25)
else:
    s=zeros(420)
    s[10]=1

ax2.plot(arange(len(s)),s,'.-')

y0=[]
for i in range(len(s)*u//d-4):
    j=i/u*d
    y0.append(s[int(j):int(j)+4]@far(j%1))
    
ax2.plot(arange(len(y0))*d/u+1,real(y0),'.-');

    
y1=signal.upfirdn(fir,s,u,d)


S=tile(fft(s),u)
FIR=fft(fir,len(S))
S*=FIR
#FIR/=max(abs(FIR))
#plot(abs(S))
#plot(abs(FIR))
#plot(abs(S*FIR))
SS=zeros(len(S)//d,'complex')
fird=fft(fir,d*128)
FIRD=zeros(len(fird)//d,'complex')
for i in range(d):
    SS+=S[i*len(S)//d:(i+1)*len(S)//d]
    FIRD+=fird[i*len(fird)//d:(i+1)*len(fird)//d]

ax1.plot(abs(fftshift(FIRD)))



y2=ifft(SS)/d
#plot(20*log10(abs(fft(fir,1024))+.01))


ax2.plot(arange(len(y1))*d/u,y1,'.-');
ax2.plot(arange(len(y2))*d/u,real(y2),'.-');
ax2.legend(['original','straigh','upfirdn','FFT'])