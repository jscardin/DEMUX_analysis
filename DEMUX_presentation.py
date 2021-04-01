from pylab import *
from matplotlib import animation
from matplotlib.widgets import Slider, Button, RadioButtons
from quicktions import Fraction

os=8   # oversampling for the entire simulation
OSfin=4  #final oversampling
span=7 # number of FFT frame, must be odd number
trim=0
adj_ch_boost=0  # adjacent channel boost

def rrc(N0,N01,N1):
    return(concatenate([zeros(N0),
                        cos(linspace(-pi,0,N01,endpoint=False))/2+0.5,
                        ones(N1),
                        cos(linspace(0,pi,N01,endpoint=False))/2+0.5,
                        zeros(N0)]))

def RRC1(N,FBW,alpha):
    '''returns zero phase RRC impulse response
    '''
    n=int(round(N/2*FBW))
    z=int(round(n*alpha))
    H=concatenate([ones(n-z),sqrt((cos(linspace(0,pi,2*z,endpoint=False))+1)/2)])
    H=concatenate([H,zeros(N-2*len(H)+1),flipud(H[1:])])
    return(H)



if True:
    fig1=figure(num=1,clear=True,facecolor='orange')

    Nfft_h         =Slider(axes([0.01,0.35,0.02,0.6]),'FFT',5,18,10,'%d',valstep=1,orientation='vertical')
    Nifft_h        =Slider(axes([0.04,0.35,0.02,0.6]),'IFFT',5,18,5,'%d',valstep=1,orientation='vertical')
    OS_h           =Slider(axes([0.07,0.35,0.02,0.6]),'OS',1,OSfin,2,'%0.2f',valstep=0.1,orientation='vertical')
    Fc_h           =Slider(axes([0.10,0.35,0.02,0.6]),'Freq',-100,100,0,'%0.1f',valstep=0.2,orientation='vertical')
    Toff_h         =Slider(axes([0.13,0.35,0.02,0.6]),'Time',-200,200,0,'%d',orientation='vertical',valstep=1)
    RollOff_h=RadioButtons(axes([0.01,0.21,.07,0.1]),('0%','5%','20%','35%'),active=1)
    OL_h     =RadioButtons(axes([0.08,0.21,.07,0.1]),('1/2','1/4','1/8','0'),active=1)
    Wtype_h  =RadioButtons(axes([0.01,0.01,.10,0.2]),('RECT OLdisc','BART OLadd','HANN OLadd','HANN OLdisc'),active=2)
    mode_h   =RadioButtons(axes([0.11,0.11,.05,0.1]),('RC-RRC','RC-FIR','POLY'),active=0)
    #RCtype_h =RadioButtons(axes([0.11,0.01,.05,0.1]),('Linear','Farrow','FLAT'),active=1)
    
    ax1=subplot2grid((8,10),(0,1),rowspan=4,colspan=9,title='Time Domain Samples',fig=fig1)
    L1a,=ax1.plot([],[],'y.-',lw=1,ms=3,label='Expected')
    L1b,=ax1.plot([],[],'r.' ,lw=1,ms=3,label='Actual')
    L1c,=ax1.plot([],[],'b.-',lw=1,ms=3,label='EVM $ \dfrac{1}{100} dB $ ')
    L1d=[ax1.plot([],[],'c--',lw=1,label='FFT Windows')[0]]
    
    for i in range(span-1):
        L1d.append(ax1.plot([],[],'c--',lw=1)[0])
    ax1.legend(loc='upper right',framealpha=1,fancybox=True)
    ax1.grid()
    ax1.set_xlabel('IFFT TIME SAMPLES')
    
    ax2=subplot2grid((8,10),(5,1),rowspan=3,colspan=9,title='Freq Domain Tones',fig=fig1)
    L2a,=ax2.plot([],[],'y.'  ,lw=2,ms=10,label='Expected')
    L2b,=ax2.plot([],[],'r.'  ,lw=1,ms=5,label='Actual')
    L2c,=ax2.plot([],[],'b'  ,lw=2,label='EVM $ \dfrac{1}{100} dB $ ')
    L2d,=ax2.plot([],[],'c--',label='Channel Envelopes')
    L2e,=ax2.plot([],[],'g',label='Interpolator Freq resp')
    ax2.legend(loc='upper right',framealpha=1,fancybox=True)
    ax2.grid()
    ax2.set_xlabel('FREQ TONES')
    
    textMERavg=fig1.text(0.3,0.04,'asdf',size='30',ha='center')
    textMER   =fig1.text(0.7,0.04,'asdf',size='30',ha='center')
    

class waveform():
    def __init__(self):
        self.calc(0)
    def calculate(self,Fc,Nfft,Nifft,OS,OL,Wtype,RollOff,mode,trim=0,plotEn=True):
        """ Extracting variables 
        """
        Nfft=int(2**Nfft)
        Nifft=int(2**Nifft)
        Toff=0
        Nbw=Nifft/OS
        OLwin=int(Nfft*eval(OL))
        OLfft=int(Nfft*eval(OL))
        OLifft=int(Nifft*eval(OL))
        if 'OLdisc' in Wtype: OLwin//=2
        ro=[a[0] for a in ((0,"0%"),(.05,"5%"),(.2,"20%"),(.35,"35%")) if a[1]==RollOff][0]
        N=Nfft*os
        """ Creating Frequency envelope
        """
        n=int(round(Nbw*os/2*(1-ro)/(1+ro)))
        Fir=sqrt(rrc(0,(int(round(Nbw*os/2))-n),2*n))
        Fenv=concatenate([Fir*10**(adj_ch_boost/20),Fir*10**(adj_ch_boost/20),Fir,Fir*10**(adj_ch_boost/20),Fir*10**(adj_ch_boost/20),[0]])
        fe=arange(len(Fenv))/os-(len(Fenv)-1)/2/os
        """ Creating main signal (s)
        """
        t=arange(-N//2,N//2)
        s=exp(1j*2*pi*Fc/Nfft*(t+Toff))*interp(Fc,fe,Fenv)
        """ Creating windows w: pre-FFT window wo: time envelope window 
        """
        if   'RECT' in Wtype:
            w=concatenate([zeros((N-Nfft)//2),ones(Nfft),zeros((N-Nfft)//2)])
            wo=concatenate([zeros((N-Nfft+OLfft)//2),ones(Nfft-OLfft),zeros((N-Nfft+OLfft)//2)])
        elif 'BART' in Wtype: 
            w=concatenate([zeros((N-Nfft)//2),
                           linspace(0,1,OLwin,endpoint=False),
                           ones(Nfft-2*OLwin),
                           linspace(1,0,OLwin,endpoint=False),
                           zeros((N-Nfft)//2)])
            wo=w
        elif 'HANN' in Wtype: 
            w=rrc((N-Nfft)//2,OLwin,Nfft-2*OLwin)
            wo=w
        """ Create overlapping FFT's s3
        """
        self.s3=zeros(span*Nifft-(span-1)*OLifft,'complex')
        f=arange(min(Nfft,Nifft*2)*os)-min(Nfft,Nifft*2)*os//2
        self.ww=zeros(N)
        for i in range(span):
            offset=(i-span//2)*(Nfft-OLfft)
            if plotEn:  L1d[i].set_data((t+offset)*Nifft/Nfft,wo[t+N//2])
            S=fftshift(fft(fftshift(w*roll(s,-offset))))/(Nfft-OLwin)
            if mode=='RC-RRC':  # select center BW only, fill with zeros outside
                n=int(ceil(Nbw/2)+trim)*2
                SS=S[arange(-n//2,n//2)*os+N//2]
            else:
                n=int(round((Nifft-Nbw)/2))
                SS=S[arange(-Nifft//2,Nifft//2)*os+N//2]*(rrc(0,n,Nifft-2*n))
            if i==span//2 and plotEn: # plot Frequency domain if center FFT
                L2a.set_data(f/os,real(S[f+N//2]))
                L2b.set_data(arange(-len(SS)//2,len(SS)//2),real(SS))
            ss=concatenate([zeros(Nifft//2-len(SS)//2),SS,zeros(Nifft//2-len(SS)//2)])
            ss=fftshift(ifft(fftshift(ss)))
            """ overlap and discard or add
            """
            if 'OLdisc' in Wtype:
                self.s3[arange(Nifft-OLifft)+i*(Nifft-OLifft)+OLifft//2]+=ss[arange(-(Nifft-OLifft)//2,(Nifft-OLifft)//2)+Nifft//2]*(Nifft-OLifft/2)
                self.ww[arange(-(Nfft-OLfft)//2,(Nfft-OLfft)//2)+offset+N//2]=1
            else:
                self.s3[arange(Nifft)+i*(Nifft-OLifft)]+=ss*(Nifft-OLifft)
                self.ww+=roll(wo,offset)
        """ Filter and interpolator
        """
        self.s3=fft(self.s3)
        if   mode=="RC-RRC": #filter first then oversample
            self.s4=self.RateConverter(self.s3)
            self.s4=self.s4*RRC1(len(self.s4),1/OSfin/(1+ro),ro)
        if   mode=="RC-FIR": #filter first then oversample
            self.s4=self.RateConverter(self.s3)
            f2=interp(linspace(-1,1,len(self.s4)),linspace(-self.up,self.up,len(self.RCfir)),fftshift(self.up/(self.RCfir+1e-20)))
            self.s4=self.s4*RRC1(len(self.s4),1/OSfin/(1+ro),ro)*fftshift(f2)
        elif mode=="POLY": #oversample first then filter
            #f2=interp(linspace(-1,1,len(self.s3)),linspace(-self.up,self.up,len(self.RCfir)),fftshift(self.up/self.RCfir))
            #self.s3=self.s3*RRC1(len(self.s3),Nbw/Nifft/(1+ro),ro)#*fftshift(f2)
            self.s4=self.RateConverter(self.s3)
        self.s4=ifft(self.s4)
        """ Creating reference signal (sr4)
        """
        sr=s*(abs(Fc)<Nbw/2)*self.ww*interp(Fc,fe,Fenv)
        i=(span*Nfft- (span-1)*OLfft)
        self.sr4=interp(linspace(-i/2,i/2,len(self.s4),endpoint=False),arange(-Nfft*os//2,Nfft*os//2),sr)
        """ Calculating MER
        """
        n=int(Nifft*OSfin/(Nifft/Nbw)*1)
        r=arange(-n,n)+len(self.s4)//2
        err=abs(self.s4[r]-self.sr4[r])**2
        """ Plotting
        """
        t4=arange(0,len(self.s3),(Nifft/Nbw)/OSfin)-len(self.s3)//2
        p4=mean(abs(self.sr4[r])**2)
        if plotEn:
            #N=16384
            #G=(1+OLfft/(span*(Nfft-OLfft)))/len(wf.s4)
            #L2a.set_data(linspace(-Nbw*OSfin/2,Nbw*OSfin/2,N,endpoint=False),fftshift(abs(fft(wf.sr4,N))*G))
            #L2b.set_data(linspace(-Nbw*OSfin/2,Nbw*OSfin/2,N,endpoint=False),fftshift(abs(fft(wf.s4 ,N))*G))
            #L2a.set_data([Fc,Fc],[sqrt(mean(abs(self.sr4[r])**2)),interp(Fc,fe,Fenv)])
            #L2b.set_data(Fc,sqrt(mean(abs(self.s4[r])**2)))
            L1a.set_data(t4,real(self.sr4))
            L1b.set_data(t4,real(self.s4))
            f=arange(-len(self.s3)//2,len(self.s3)//2)
            L2e.set_data(f/len(self.s3)*Nifft,fftshift(self.RCfir)[f+len(self.RCfir)//2]/self.up)
            ax1.axis([-Nifft*os//2,Nifft*os//2,-1.5,1.5])
            textMER.set_text('MER spot=%0.1f'%(-10*log10(mean(err)+1e-20)))
            L2d.set_data(fe,Fenv)
            ax2.set_xticks([-Nifft//2,-Nbw/2,Nbw/2,Nifft//2])
            ax2.axis([-Nifft*0.6,Nifft*0.6,-0.5,1.5])
            fig1.canvas.draw()
            #fig1.canvas.flush_events()
        return(t4[r],err,p4)
    """_______________________________________________________________
    """

    def RateConverter(self,s):
        s=tile(s,self.up)*self.RCfir
        z=int(ceil(len(s)/self.dn))*self.dn-len(s)
        ss=concatenate([s[:len(s)//2],zeros(z),s[-len(s)//2:]])
        y1=sum(reshape(ss,[self.dn,-1]),0)/self.dn
        return(y1)
    def SweepFreq(self,Nfft,Nifft,OS,OL,Wtype,RollOff,mode,trim=0,PlotEn=True):
        MERfreq=[]
        N=200
        Fsweep=linspace(-2**Nifft,2**Nifft,N)
        mm=array([])
        pp=0
        for Fc in Fsweep:
            t,m,p=self.calculate(Fc,Nfft,Nifft,OS,OL,Wtype,RollOff,mode,trim,False)
            MERfreq.append(-10*log10(mean(m)+1e-20)/100)
            mm=mm+m if len(mm) else m
            pp+=p
            #print('%d %0.1f %0.2f'%(Fc,10*log10(mean(m)),pp))
        self.MERavg=-10*log10(mean(mm)/pp)
        if PlotEn:
            L1c.set_data(t,-10*log10(mm/pp+1e-20)/100)
            L2c.set_data(Fsweep,MERfreq)
            textMERavg.set_text('MER avg=%0.1fdB'%(self.MERavg))
        return(self.MERavg)
    def IFFT_calc(self,Nfft,Nifft,OS,OL,Wtype,RollOff,mode,trim=0,PlotEn=True):
        ro=eval(RollOff[:-1])/100
        OLifft=int(2**Nifft*eval(OL))
        n=int(span*2**Nifft-(span-1)*OLifft)
        (self.up,self.dn)=Fraction(OSfin/OS).limit_denominator().as_integer_ratio()
        if mode=="POLY":
            #self.RCfir=fftshift(concatenate([zeros(n*self.up//2-n//2-1),ones(n+2),zeros(n*self.up//2-n//2-1)]))*self.up     
            self.RCfir=(RRC1(n*self.up,1/self.up/(1+ro)/OS,ro))*self.up
        else:
            #far=lambda mu:array([1-mu,mu])   # Linear interpolation
            far=lambda mu:array([mu*(-0.5+mu/2),-mu/2-mu**2/2+1,mu*(1.5-mu/2),-mu/2+mu**2/2]) #Farrow interpolation
            fir=reshape(list(zip(*[far(i/self.up) for i in range(self.up,0,-1)])),[-1])
            self.RCfir=real(fft(fftshift(concatenate([zeros(n*self.up//2-len(fir)//2),fir,zeros(n*self.up//2-len(fir)//2)]))))  
        mer=self.SweepFreq(   Nfft,Nifft,OS,OL,Wtype,RollOff,mode,trim,PlotEn)
        print('config: ',Nfft,Nifft,'%0.2f'%OS,OL,Wtype,RollOff,mode,mer)
        return(mer)
    def calc(self,a):
        self.IFFT_calc(         Nfft_h.val,Nifft_h.val,OS_h.val,OL_h.value_selected,Wtype_h.value_selected,RollOff_h.value_selected,mode_h.value_selected,trim)
        self.calculate(Fc_h.val,Nfft_h.val,Nifft_h.val,OS_h.val,OL_h.value_selected,Wtype_h.value_selected,RollOff_h.value_selected,mode_h.value_selected,trim)
    def calc_fc(self,Fc):
        self.calculate(Fc,      Nfft_h.val,Nifft_h.val,OS_h.val,OL_h.value_selected,Wtype_h.value_selected,RollOff_h.value_selected,mode_h.value_selected,trim)
            
wf=waveform()

Fc_h.on_changed(wf.calc_fc)

Toff_h.on_changed(wf.calc)
Nfft_h.on_changed(wf.calc)
Nifft_h.on_changed(wf.calc)
OS_h.on_changed(wf.calc)

RollOff_h.on_clicked(wf.calc)
OL_h.on_clicked(wf.calc)
Wtype_h.on_clicked(wf.calc)
mode_h.on_clicked(wf.calc)
#RCtype_h.on_clicked(wf.calc)

show()


