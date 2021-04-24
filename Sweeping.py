from pylab import *
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time



fig2,ax22=subplots(1,1,num=2,clear=True)

#Nfft=2**13;OL=0.5;Wtype='RECT OLdisc';alpha=.05;mode='POLY';FIRlen=384*2;PolyWidth=4;OSrange=arange(4,2,-0.05);RRCbeta=0
Nfft=2**13;OL=0.25;Wtype='RECT OLdisc';alpha=.05;mode='RC-RRC';FIRlen=31;OSrange=arange(3,1.5,-0.05);RRCbeta=0;PolyWidth=0
#Nfft=2**13;OL=0.5;Wtype='RECT OLdisc';alpha=.05;mode='RC-FIR';FIRlen=[31,65,129,257,513];OSrange=arange(4,2,-0.2);RRCbeta=2
#Nfft=2**13;OL=0.5;Wtype='RECT OLdisc';alpha=.05;mode='RC-FIR';FIRlen=128;OSrange=arange(4,2,-0.02);RRCbeta=2


for i in range(len(fig2.texts)):
    fig2.texts[0].remove()

y=0.97;yy=0.04;x=0;xx=0.15
fig2.text(x,y,'Window Type= '+Wtype);y-=yy
fig2.text(x,y,'OS range %g to %g'%(min(OSrange),max(OSrange)));y-=yy
fig2.text(x,y,'OL ratio=1/%d'%(1/OL));y-=yy
fig2.text(x,y,'Roll Off=%d%%'%(alpha*100));y-=yy
fig2.text(x,y,'RRCbeta=%g'%RRCbeta);y-=yy
#fig2.text(x,y,'mode: %s'%mode);y-=yy
#fig2.text(x,y,'Nfft=%d'%Nfft);y-=yy
#fig2.text(10,y,'FIR len=%d'%(FIRlen));y-=yy
#fig2.text(10,y,'PolyWidth=%d'%(PolyWidth));y-=yy

ax22.grid()
ax22.set_ylabel('EVM (dB)')
ax22.set_xlabel('Number of Tones per Channel')
#ax22.set_xscale('log')

t1=time.time()
colors=['r','g','b','c','y']*10
leg=[]

#for mode in ['RC-RRC','RC-FIR','POLY','PERFECT']:
#for (mode,FIRlen) in [('RC-RRC',31),('RC-RRC',61),('RC-FIR',61),]:
#for Wtype in ['RECT OLdisc','HANN OLadd']:
#for OL in [0.125]:
#for PolyWidth in [2,4,6]:
#for RRCbeta in [0,1,2,3,4]:
#for FIRlen in FIRlen:
#for (Wtype,OL) in [('HANN OLadd',0.25),('RECT OLdisc',0.5)]:
for a in [1]:
    #co=colors.pop()
    #ty=['.-','.--']*10
#    if mode=='RC-RRC':
#        OSrange=arange(3,1.5,-0.05)
#        FIRlen=31
#        Wtype='RECT OLdisc'
    tones=[]
    mer=[]
    #leg.append('RRCbeta=%d  Nifft=%d'%(RRCbeta,Nifft))
    leg.append('MODE=%s FIRlen=%d'%(mode,FIRlen))
    #leg.append('OL type=1/%d %s'%(1/OL,Wtype))
    print('_____',leg[-1])
    with ProcessPoolExecutor(max_workers = 16) as e:
        for Nifft in [32,64,128]:
            for OS in OSrange:
                mer.append(e.submit(wf.standalone,Nfft,Nifft,OS,OL,alpha,Wtype,mode,FIRlen,PolyWidth,RRCbeta))
            tones=concatenate([tones,Nifft/OSrange])
    ax22.plot(tones,[i.result() for i in mer])
    fig2.canvas.draw()
    fig2.canvas.flush_events()
        

ax22.legend(leg)
print('elapse time (mins)',(time.time()-t1)/60)