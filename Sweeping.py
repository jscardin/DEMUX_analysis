from pylab import *
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time



fig2,ax22=subplots(1,1,num=2,clear=True)

Nfft=2**13
Nifft=2**5
OL=0.5
Wtype='HANN OLadd'
alpha=.05
mode='RC-FIR'
FIRlen=511
OSrange=arange(4,2,-0.05)


y=80;yy=5
ax22.axis([8,850,0,y+10])
ax22.text(10,y,'Nfft=%d'%Nfft);y-=yy
ax22.text(10,y,'OL ratio=1/%d'%(1/OL));y-=yy
#ax22.text(10,y,'mode: %s'%mode);y-=yy
ax22.text(10,y,'Roll Off=%d%%'%(alpha*100));y-=yy
ax22.text(10,y,'FIR len=%d'%(FIRlen));y-=yy
ax22.grid()
ax22.set_ylabel('EVM (dB)')
ax22.set_xlabel('Number of Tones per Channel')
ax22.set_xscale('log')

t1=time.time()
colors=['r','g','b','c','y']*10
leg=[]

#for mode in ['RC-FIR','POLY']:
for FIRlen in [31,63,127,255,511]:
    co=colors.pop()
    ty=['.-','.--']*10
#    if mode=='RC-RRC':
#        OSrange=arange(3,1.5,-0.05)
#        FIRlen=31
#        Wtype='RECT OLdisc'
    
    for Nifft in [32,64,128,256,512]:
        leg.append('FIRlen=%s  Nifft=%d'%(FIRlen,Nifft))
        print('_____',leg[-1])
    
        mer=[]
        with ProcessPoolExecutor(max_workers = 16) as e:
            for OS in OSrange:
                mer.append(e.submit(wf.standalone,Nfft,Nifft,OS,OL,alpha,Wtype,mode,FIRlen))
        mer=[i.result() for i in mer]
        ax22.plot(Nifft/OSrange,mer,co+ty.pop())
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        

ax22.legend(leg)
print('elapse time ',time.time()-t1)