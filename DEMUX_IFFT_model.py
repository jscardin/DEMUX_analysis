

from pylab import *

fig=figure(num=1,clear=True)

Nifft_h        =Slider(axes([0.04,0.35,0.02,0.6]),'IFFT',5,16,7,'%d',valstep=1,orientation='vertical')
OS_h           =Slider(axes([0.07,0.35,0.02,0.6]),'OS',1,4,2,'%0.2f',valstep=0.05,orientation='vertical')
#Fc_h           =Slider(axes([0.10,0.35,0.02,0.6]),'Freq',-100,100,2,'%0.1f',valstep=0.2,orientation='vertical')
Toff_h         =Slider(axes([0.13,0.35,0.02,0.6]),'Time',-1,1,0,'%d',orientation='vertical')
# RollOff_h=RadioButtons(axes([0.01,0.21,.07,0.1]),('0%','5%','20%','35%'),active=1)
# OL_h     =RadioButtons(axes([0.08,0.21,.07,0.1]),('1/2','1/4','1/8','0'),active=1)
# Wtype_h  =RadioButtons(axes([0.01,0.01,.10,0.2]),('RECT OL discard','BARTLETT OL add','HANN OL add','HANN OL discard'),active=2)
# mode_h   =RadioButtons(axes([0.11,0.01,.05,0.1]),('IDFT','FIR'),active=0)



axt=subplot2grid((8,10),(0,1),rowspan=4,colspan=9,title='DEMUX IMPULSE RESPONSE')
axf=subplot2grid((8,10),(5,1),rowspan=3,colspan=9,title='DEMUX FREQ RESPONSE')

L1finite,=axt.plot([],[],'.-',lw=1,label='Impulse response finite')
L1infini,=axt.plot([],[],'.-',lw=1,label='Impulse response infinite')
L2finite,=axf.plot([],[],lw=2,label='Freq response finite')
L2infini,=axf.plot([],[],'.-',lw=2,label='Freq response infinite')

axf.grid()
axt.grid()


class waveform():
    def __init__(self):
        self.calculate(1)

    def calculate(self,a):
        Nifft=int(2**Nifft_h.val)
        N=Nifft*16
        BW=int(Nifft/OS_h.val)
        R=int(Toff_h.val*Nifft/2)
        
        H=zeros(Nifft//2+1)
        
        H[:BW//2]=1
        h=fftshift(irfft(H))
        h=roll(h,R)
        H=fft(h)
        hh=concatenate([zeros(N//2-Nifft//2),h,zeros(N//2-Nifft//2)])
        
        HH=fft(hh)
        
        
        L1finite.set_data(arange(-N//2,N//2),hh*Nifft/BW)
        L1infini.set_data(arange(-Nifft//2,Nifft//2),h*Nifft/BW)
        
        L2finite.set_data(linspace(-Nifft//2,Nifft//2,N),20*log10(abs(fftshift(HH))+1e-20))
        L2infini.set_data(arange(-Nifft//2,Nifft//2),20*log10(abs(fftshift(H))+1e-20))
        
        axt.axis([-Nifft//2-10,Nifft//2+10,-1,1])
        axf.axis([-Nifft//2,Nifft//2,-40,5])
        fig.canvas.draw()
        #fig.canvas.flush_events()
    def rrc(N0,N01,N1):
        return(concatenate([zeros(N0),
                        cos(linspace(-pi,0,N01,endpoint=False))/2+0.5,
                        ones(N1),
                        cos(linspace(0,pi,N01,endpoint=False))/2+0.5,
                        zeros(N0)]))


wf=waveform()

#Nfft_h.on_changed(wf.calculate)
Nifft_h.on_changed(wf.calculate)
#Fc_h.on_changed(wf.calculate)
Toff_h.on_changed(wf.calculate)
OS_h.on_changed(wf.calculate)
#RollOff_h.on_clicked(wf.calculate)
#OL_h.on_clicked(wf.calculate)
#Wtype_h.on_clicked(wf.calculate)
#mode_h.on_clicked(wf.calculate)

