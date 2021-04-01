from pylab import *
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time



fig2,ax22=subplots(1,1,num=2,clear=True)

Nfft=13
Nifft=5
OL='1/4'
Wtype='HANN OLadd'
RollOff='5%'
mode='RC-RRC'


y=95;yy=5
ax22.text(10,y,'Nfft=%d'%(2**Nfft));y-=yy
ax22.text(10,y,'OL ratio=%s'%(OL));y-=yy
ax22.text(10,y,'OL mode: %s'%Wtype);y-=yy
ax22.text(10,y,'Roll Off=%s'%RollOff);y-=yy
#fig2.text(0,y,'Filter mode=%s'%mode);y-=yy
ax22.grid()
ax22.set_ylabel('EVM (dB)')
ax22.set_xlabel('Number of Tones per Channel')
ax22.set_xscale('log')
ax22.axis([8,512,0,100])

t1=time.time()
colors=['r','g','b','c','y']*10
leg=[]

for mode in [a.get_text() for a in mode_h.labels]:
    co=colors.pop()
    ty=['.-','.--']*10
    for Nifft in [5,6,7,8,9]:
        leg.append('%s  Nifft=%d'%(mode,2**Nifft))
        print('_____',leg[-1])
        if mode=='RC-RRC':
            OSrange=arange(3,1.5,-0.05)
        else:
            OSrange=arange(4,2,-0.05)
        mer=[]
        with ProcessPoolExecutor(max_workers = 16) as e:
            for OS in OSrange:
                mer.append(e.submit(wf.IFFT_calc,Nfft,Nifft,OS,OL,Wtype,RollOff,mode,False))
        mer=[i.result() for i in mer]
        ax22.plot(2**Nifft/OSrange,mer,co+ty.pop())
        fig2.canvas.draw()
        fig2.canvas.flush_events()
    

ax22.legend(leg)
print('elapse time ',time.time()-t1)