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
    H=concatenate([ones(n-z),(sqrt(cos(linspace(0,pi,2*z,endpoint=False))+1)/2),zeros(NN//2+1-n-z)])
    h=irfft(H)
    h=concatenate([h[-N//2:],h[:N//2]])#*kaiser(N,beta)
    return(h)

far=lambda mu:array([1-mu,mu])
far=lambda mu:array([mu*(-0.5+mu/2),-mu/2-mu**2/2+1,mu*(1.5-mu/2),-mu/2+mu**2/2])

def interpFIR(u,d,n):
    fir=reshape(list(zip(*[far(i/u) for i in range(u,0,-1)])),[-1])
    FIR=fft(fftshift(concatenate([zeros(n*u//2-len(fir)//2),fir,zeros(n*u//2-len(fir)//2)])))
    FIRds=real(concatenate([FIR[-n*u//d//2:],FIR[:n*u//d//2]]))/u
    f=linspace(-u*n/d/2,u*n/d/2,len(FIRds))
    return(f,FIRds)


#fig2,ax=subplots(1,num=2,clear=True)

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
L1a,=ax1.plot([],[],'y.-',lw=1,label='real(Input Signal)')
L1b,=ax1.plot([],[],'b.',lw=1,label='real(DEMUX output)')
L1c=[ax1.plot([],[],'r',lw=1,label='FFT blocks')[0]]
for i in range(span-1):
    L1c.append(ax1.plot([],[],'r',lw=1)[0])
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_xlabel('IFFT TIME SAMPLES')


ax2=subplot2grid((8,10),(5,1),rowspan=3,colspan=9,title='Freq Domain Tones')
L2a,=ax2.plot([],[],'y',lw=4,label='real(Input Signal)')
L2b,=ax2.plot([],[],'b.',label='real(IFFT input)')
L2c,=ax2.plot([],[],'r',label='Channel Envelope')
L2d,=ax2.plot([],[],label='Interpolator Freq resp')
ax2.legend()
ax2.grid()
ax2.set_xlabel('IFFT FREQ TONES')



textMER=ax1.text(0,-1.5,'asdf')


class waveform():
    def __init__(self):
        self.plotEn=True

    def SweepFreq(self,a):
        self.plotEn=False
        for f in arange(-100,100):
            print(self.calculate(f))
        self.plotEn=True


    def calculate(self,Fcs):
        #______ Extracting variables
        Nfft=int(2**Nfft_h.val)
        Nifft=int(2**Nifft_h.val)
        Nbw=Nifft/OS_h.val
        OLwin=int(Nfft*eval(OL_h.value_selected))
        OLfft=int(Nfft*eval(OL_h.value_selected))
        OLifft=int(Nifft*eval(OL_h.value_selected))
        if 'discard' in Wtype_h.value_selected: OLwin//=2
        Fc=Fcs*Nifft/160
        ro=[a[0] for a in ((0,"0%"),(.05,"5%"),(.2,"20%"),(.35,"35%")) if a[1]==RollOff_h.value_selected][0]
        N=Nfft*os
        #______ Creating Frequency envelope
        Fir=rrc(0,int(round(os*Nbw*ro/2)),int(round(os*Nbw*(1-ro))))
        Fenv=concatenate([tile(Fir,span),[0]])
        fe=arange(len(Fenv))/os-(len(Fenv)-1)/2/os
        #______ Creating main signal (s)
        t=arange(-N//2,N//2)
        s=exp(1j*2*pi*Fc/Nfft*(t+Toff_h.val))*interp(Fc,fe,Fenv)
        #______ Creating windows w: pre-FFT window wo: time envelope window 
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
        #______ Create overlapping FFT's s3        
        self.s3=zeros(span*Nifft-(span-1)*OLifft,'complex')
        f=arange(min(Nfft,Nifft*2)*os)-min(Nfft,Nifft*2)*os//2
        self.ww=zeros(N)
        for i in range(span):
            offset=(i-span//2)*(Nfft-OLfft)
            if self.plotEn:  L1c[i].set_data((t+offset)*Nifft/Nfft,wo[t+N//2])
            S=fftshift(fft(fftshift(w*roll(s,-offset))))/(Nfft-OLwin)
            if mode_h.value_selected=="FIR": #Select whole range
                SS=S[arange(-Nifft//2,Nifft//2)*os+N//2]
            else: # select center BW only, fill with zeros outside
                n=int(round(Nbw/2)*2)
                SS=S[arange(-n//2,n//2)*os+N//2]
            if i==span//2 and self.plotEn: # plot Frequency domain if center FFT
                L2a.set_data(f/os,real(S[f+N//2]))
                L2b.set_data(arange(-len(SS)//2,len(SS)//2),real(SS))
            ss=concatenate([zeros(Nifft//2-len(SS)//2),SS,zeros(Nifft//2-len(SS)//2)])
            ss=fftshift(ifft(fftshift(ss)))
            #______ overlap and discard or add
            if 'discard' in Wtype_h.value_selected:
                self.s3[arange(Nifft-OLifft)+i*(Nifft-OLifft)+OLifft//2]+=ss[arange(-(Nifft-OLifft)//2,(Nifft-OLifft)//2)+Nifft//2]*(Nifft-OLifft/2)
                self.ww[arange(-(Nfft-OLfft)//2,(Nfft-OLfft)//2)+offset+N//2]=1
            else:
                self.s3[arange(Nifft)+i*(Nifft-OLifft)]+=ss*(Nifft-OLifft)
                self.ww+=roll(wo,offset)
 
        #______ Filter and interpolator
        t4=arange(0,len(self.s3),(Nifft/Nbw)/OSfin)    
        t4=t4[:len(t4)//2*2]
        if mode_h.value_selected=="FIR": #filter first then oversample
            self.s3=self.RRCfilter(self.s3,Nbw/Nifft,ro)
            self.s4=self.RateConverter(t4,self.s3)
        else: #oversample first then filter
            self.s4=self.RateConverter(t4,self.s3)
            self.s4=self.RRCfilter(self.s4,1/OSfin,ro)
        #______ Creating reference signal (sr4)
        sr=s*(abs(Fc)<Nbw/2)*self.ww
        i=(span*Nfft- (span-1)*OLfft)
        j=(span*Nifft-(span-1)*OLifft)
        self.sr4=interp(linspace(-i//2,i//2,len(self.s4),endpoint=False),arange(-Nfft*os//2,Nfft*os//2),sr)
        #______ Calculating MER
        MER=10*log10(sum(abs(self.s4-self.sr4)**2)+1e-20)
        #______ Plotting
        if self.plotEn:
            t4=t4-len(self.s3)//2
            L1b.set_data(t4,real(self.s4))
            L2d.set_data(interpFIR(int(OSfin*Nbw),Nifft,Nifft))
            L1a.set_data(linspace(-j//2,j//2,len(self.sr4),endpoint=False),real(self.sr4))
            ax1.axis([-Nifft*os//2,Nifft*os//2,-1.5,1.5])
            textMER.set_text('MER est.=%0.1f'%(MER))
            L2c.set_data(fe,Fenv)
            ax2.set_xticks([-Nifft//2,-Nbw/2,Nbw/2,Nifft//2])
            ax2.axis([-Nifft*0.6,Nifft*0.6,-0.5,1.5])
            fig.canvas.draw()
            #fig.canvas.flush_events()
        return(MER)

    def RateConverter(self,t,s):
        y0=[]
        lenfar=len(far(0))
        for j in roll(t,0):
            y0.append(s[arange(int(j)-1,int(j)+3)%len(s)]@far(j%1))
        return(roll(array(y0),0))
#        return(interp(t,arange(len(s)),s))
    def RRCfilter(self,s,bw,ro):
            n=len(s)//2
            ss=roll(convolve(s,RRC(FIRlen,bw,ro,2),'same'),-1)
            return(ss[-n+len(ss)//2:n+len(ss)//2]) #extract center part only if FIR too long
    def button_pressed(self,Nifft):
        ax1.axis([-2**Nifft*os//2,2**Nifft*os//2,-1.5,1.5])
        fig.canvas.toolbar.update()
        self.calculate(0)

wf=waveform()


#fig.canvas.mpl_connect('button_press_event',wf.button_pressed)

Fc_h.on_changed(wf.calculate)


Nfft_h.on_changed(wf.calculate)
Nifft_h.on_changed(wf.calculate)
Toff_h.on_changed(wf.calculate)
OS_h.on_changed(wf.calculate)
RollOff_h.on_clicked(wf.calculate)
OL_h.on_clicked(wf.calculate)
Wtype_h.on_clicked(wf.calculate)
mode_h.on_clicked(wf.calculate)

show()