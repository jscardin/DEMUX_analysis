from pylab import *
from matplotlib import animation
from matplotlib.widgets import Slider, Button, RadioButtons

os=8   # oversampling for the entire simulation
OSfin=4  #final oversampling
span=7 # number of FFT frame, must be odd number

def rrc(N0,N01,N1):
    return(concatenate([zeros(N0),
                        cos(linspace(-pi,0,N01,endpoint=False))/2+0.5,
                        ones(N1),
                        cos(linspace(0,pi,N01,endpoint=False))/2+0.5,
                        zeros(N0)]))

def RRC(N,FBW,alpha,beta=0):
    ''' returns zero phase RRC impulse response'''
    NN=N*8
    n=int(round(NN/2*FBW))
    z=int(round(n/2*alpha))
    H=concatenate([ones(n-z),sqrt((cos(linspace(0,pi,2*z,endpoint=False))+1)/2),zeros(NN//2+1-n-z)])
    h=irfft(H)
    h=concatenate([h[-N//2:],h[:N//2]])*kaiser(N,beta)
    return(h)


fig2,ax=subplots(1,num=2,clear=True)

#fig,(a,ax1,ax2)=subplots(1,3,num=1,clear=True)
fig=figure(num=1,clear=True)

Nfft_h         =Slider(axes([0.01,0.35,0.02,0.6]),'FFT',5,18,12,'%d',valstep=1,orientation='vertical')
Nifft_h        =Slider(axes([0.04,0.35,0.02,0.6]),'IFFT',5,18,7,'%d',valstep=1,orientation='vertical')
OS_h           =Slider(axes([0.07,0.35,0.02,0.6]),'OS',1,OSfin,2,'%0.2f',valstep=0.25,orientation='vertical')
Fc_h           =Slider(axes([0.10,0.35,0.02,0.6]),'Freq',-100,100,2,'%0.1f',valstep=0.2,orientation='vertical')
Toff_h         =Slider(axes([0.13,0.35,0.02,0.6]),'Time',-200,200,0,'%d',orientation='vertical',valstep=1)
RollOff_h=RadioButtons(axes([0.01,0.21,.07,0.1]),('0%','5%','20%','35%'),active=1)
OL_h     =RadioButtons(axes([0.08,0.21,.07,0.1]),('1/2','1/4','1/8','0'),active=1)
Wtype_h  =RadioButtons(axes([0.01,0.01,.10,0.2]),('RECT OL discard','BARTLETT OL add','HANN OL add','HANN OL discard'),active=2)
mode_h   =RadioButtons(axes([0.11,0.01,.05,0.1]),('IDFT','FIR'),active=0)
FIRlen=512


ax1=subplot2grid((8,10),(0,1),rowspan=4,colspan=9,title='Time Domain Samples')
L1InReal,=ax1.plot([],[],'y.-',lw=1,label='real(Input Signal)')
L1OutReal,=ax1.plot([],[],'b.',lw=1,label='real(DEMUX output)')
L1win=[ax1.plot([],[],'r',lw=1,label='FFT blocks')[0]]
for i in range(span-1):
    L1win.append(ax1.plot([],[],'r',lw=1)[0])
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_xlabel('IFFT TIME SAMPLES')


ax2=subplot2grid((8,10),(5,1),rowspan=3,colspan=9,title='Freq Domain Tones')
L2fftoutreal,=ax2.plot([],[],'y',lw=4,label='real(Input Signal)')
L2fftpoints,=ax2.plot([],[],'b.',label='real(IFFT input)')
L2env,=ax2.plot([],[],'r',label='Channel Envelope')
ax2.legend()
ax2.grid()
ax2.set_xlabel('IFFT FREQ TONES')



textMER=ax1.text(0,-1.5,'asdf')


class waveform():
    def __init__(self):
        self.calculate(1)


    def calculate(self,a):
        Nfft=int(2**Nfft_h.val)
        Nifft=int(2**Nifft_h.val)
        self.Nbw=Nifft/OS_h.val
        OLwin=int(Nfft*eval(OL_h.value_selected))
        OLfft=int(Nfft*eval(OL_h.value_selected))
        OLifft=int(Nifft*eval(OL_h.value_selected))
        if 'discard' in Wtype_h.value_selected: OLwin//=2
        Fc=Fc_h.val*Nifft/160
        ro=[a[0] for a in ((0,"0%"),(.05,"5%"),(.2,"20%"),(.35,"35%")) if a[1]==RollOff_h.value_selected][0]
        N=Nfft*os
        Fir=rrc(0,int(round(os*self.Nbw*ro/2)),int(round(os*self.Nbw*(1-ro))))
        Fenv=concatenate([tile(Fir,span),[0]])
        print('Fenv=%d'%len(Fenv))
        f=arange(len(Fenv))/os-(len(Fenv)-1)/2/os
        L2env.set_data(f,Fenv)
        t=arange(-N//2,N//2)
        s=exp(1j*2*pi*Fc/Nfft*(t+Toff_h.val))*interp(Fc,f,Fenv)
        if 'RECT' in Wtype_h.value_selected:
            w=concatenate([zeros((N-Nfft)//2),ones(Nfft),zeros((N-Nfft)//2)])
            wo=concatenate([zeros((N-Nfft+OLfft)//2),ones(Nfft-OLfft),zeros((N-Nfft+OLfft)//2)])
        elif 'BARTLETT' in Wtype_h.value_selected: 
            w=concatenate([zeros((N-Nfft)//2),
                           linspace(0,1,OLwin,endpoint=False),
                           ones(Nfft-2*OLwin),
                           linspace(1,0,OLwin,endpoint=False),
                           zeros((N-Nfft)//2)])
            wo=w
        elif 'HANN' in Wtype_h.value_selected: 
            w=rrc((N-Nfft)//2,OLwin,Nfft-2*OLwin)
            wo=w
        
        
        self.s3=zeros(span*Nifft-(span-1)*OLifft,'complex')
        f=arange(min(Nfft,Nifft*2)*os)-min(Nfft,Nifft*2)*os//2
        #L2fir.set_data(f,real(fftshift(FIR4)))
        ww=zeros(N)
        for i in range(span):
            ov=(i-span//2)*(Nfft-OLfft)
            L1win[i].set_data((t+ov)*Nifft/Nfft,wo[t+N//2])
            S=fftshift(fft(fftshift(w*roll(s,-ov))))/(Nfft-OLwin)
            if i==span//2: L2fftoutreal.set_data(f/os,real(S[f+N//2]))
            if mode_h.value_selected=="FIR":
                SS=S[arange(-Nifft//2,Nifft//2)*os+N//2]
            else:
                n=int(round(self.Nbw/2)*2)
                SS=S[arange(-n//2,n//2)*os+N//2]
            if i==span//2: L2fftpoints.set_data(arange(-len(SS)//2,len(SS)//2),real(SS))
            ss=concatenate([zeros(Nifft//2-len(SS)//2),SS,zeros(Nifft//2-len(SS)//2)])
            ss=fftshift(ifft(fftshift(ss)))

            if 'discard' in Wtype_h.value_selected:
                self.s3[arange(Nifft-OLifft)+i*(Nifft-OLifft)+OLifft//2]+=ss[arange(-(Nifft-OLifft)//2,(Nifft-OLifft)//2)+Nifft//2]*(Nifft-OLifft/2)
                ww[arange(-(Nfft-OLfft)//2,(Nfft-OLfft)//2)+ov+N//2]=1
            else:
                self.s3[arange(Nifft)+i*(Nifft-OLifft)]+=ss*(Nifft-OLifft)
                ww+=roll(wo,ov)
        
        if mode_h.value_selected=="FIR":
            self.s3=roll(convolve(self.s3,RRC(FIRlen,self.Nbw/Nifft,ro,2),'same'),-1)
        else:
            self.s4=interp(linspace(-len(self.s3)//2,len(self.s3)//2,OSfin/))
            pass
        
        sr=s*(abs(Fc)<self.Nbw/2)*ww
        L1InReal.set_data(t*Nifft/Nfft,real(sr))
        L1OutReal.set_data(arange(len(self.s3))-len(self.s3)//2,real(self.s3))
        i=arange(-Nifft//2,Nifft//2)*Nfft//Nifft+N//2
        j=arange(-Nifft//2,Nifft//2)+len(self.s3)//2
        textMER.set_text('MER est.=%0.1f'%(10*log10(sum(abs(self.s3[j]-sr[i])**2)/sum(abs(ww[i])**2)+1e-20)))
        ax1.axis([-Nifft*2.2,Nifft*2.2,-1.5,1.5])
        #a=array([-Nfft*3//2+OLifft,-Nfft*3//2+2*OLifft,-Nfft*3//2+OLifft,-Nfft//2,-Nfft//2+OLifft])
        #ax1.set_xticks(concatenate([a,-a]))
        ax2.set_xticks([-Nifft//2,-self.Nbw/2,self.Nbw/2,Nifft//2])
        ax2.axis([-Nifft*0.6,Nifft*0.6,-0.5,1.5])
        fig.canvas.draw()
        #fig.canvas.flush_events()
        

wf=waveform()

fig.canvas.toolbar.update()

Nfft_h.on_changed(wf.calculate)
Nifft_h.on_changed(wf.calculate)
Fc_h.on_changed(wf.calculate)
Toff_h.on_changed(wf.calculate)
OS_h.on_changed(wf.calculate)
RollOff_h.on_clicked(wf.calculate)
OL_h.on_clicked(wf.calculate)
Wtype_h.on_clicked(wf.calculate)
mode_h.on_clicked(wf.calculate)

