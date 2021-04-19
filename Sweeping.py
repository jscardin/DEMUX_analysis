from pylab import *
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time



fig2,ax22=subplots(1,1,num=2,clear=True)

Nfft=2**13;Nifft=2**5;OL=0.5;Wtype='RECT OLdisc';alpha=.05;mode='POLY';FIRlen=384*2;PolyWidth=4;OSrange=arange(4,2,-0.05);RRCbeta=0
Nfft=2**13;Nifft=2**5;OL=0.25;Wtype='RECT OLdisc';alpha=.05;mode='RC-RRC';FIRlen=31;OSrange=arange(3,1.5,-0.2);RRCbeta=0


y=10;yy=2
ax22.axis([8,850,0,y+10])
ax22.text(10,y,'Nfft=%d'%Nfft);y-=yy
ax22.text(10,y,'OL ratio=1/%d'%(1/OL));y-=yy
ax22.text(10,y,'Window Type= '+Wtype);y-=yy
ax22.text(10,y,'mode: %s'%mode);y-=yy
ax22.text(10,y,'Roll Off=%d%%'%(alpha*100));y-=yy
#ax22.text(10,y,'RRCbeta=%0.1f'%RRCbeta);y-=yy
#ax22.text(10,y,'FIR len=%d'%(FIRlen));y-=yy
#ax22.text(10,y,'PolyWidth=%d'%(PolyWidth));y-=yy

ax22.grid()
ax22.set_ylabel('EVM (dB)')
ax22.set_xlabel('Number of Tones per Channel')
ax22.set_xscale('log')

t1=time.time()
colors=['r','g','b','c','y']*10
leg=[]

#for mode in ['RC-RRC','RC-FIR','POLY','PERFECT']:
#for Wtype in ['RECT OLdisc','HANN OLadd']:
#for OL in [0.125]:
#for PolyWidth in [2,4,6]:
for RRCbeta in [0,1,2,3,4]:
#for FIRlen in [31,61]:
    co=colors.pop()
    ty=['.-','.--']*10
#    if mode=='RC-RRC':
#        OSrange=arange(3,1.5,-0.05)
#        FIRlen=31
#        Wtype='RECT OLdisc'
    
    for Nifft in [32,64,128,256]:
        leg.append('RRCbeta=%d  Nifft=%d'%(RRCbeta,Nifft))
        print('_____',leg[-1])
    
        mer=[]
        with ProcessPoolExecutor(max_workers = 16) as e:
            for OS in OSrange:
                mer.append(e.submit(wf.standalone,Nfft,Nifft,OS,OL,alpha,Wtype,mode,FIRlen,PolyWidth,RRCbeta))
        mer=[i.result() for i in mer]
        ax22.plot(Nifft/OSrange,mer,co+ty.pop())
        fig2.canvas.draw()
        fig2.canvas.flush_events()
        

ax22.legend(leg)
print('elapse time (mins)',(time.time()-t1)/60)