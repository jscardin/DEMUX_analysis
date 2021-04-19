from pylab import *
from matplotlib import animation
from matplotlib.widgets import Slider, Button, RadioButtons,TextBox
from quicktions import Fraction

os=12   # oversampling for the entire simulation
OSfin=4  #final oversampling
span=7 # number of FFT frame, must be odd number
trim=0
Nchan=3
adj_ch_boost=0  # adjacent channel boost
NpointFsweep=200 # number points in the frequency sweep.


def rrc(N0,N01,N1):
    return(concatenate([zeros(N0),
                        cos(linspace(-pi,0,N01,endpoint=False))/2+0.5,
                        ones(N1),
                        cos(linspace(0,pi,N01,endpoint=False))/2+0.5,
                        zeros(N0)]))

def RRC1(N,FIRlen,FBW,alpha,Beta=0):
    n=int(min(FIRlen,N))//2
    t=arange(-n,n+1)*FBW+1e-10
    Hf = FBW*(sin(pi*t*(1-alpha))+(4*alpha*t*cos(pi*t*(1+alpha))))/(pi*t*(1-(4*alpha*t)**2))
    Hf*=kaiser(len(Hf),Beta)
    hf=concatenate([Hf[-n-1:],zeros(N-len(Hf)),Hf[:n]])
    return(fft(hf/sum(hf)))

#    
#def RRC1(N,FIRlen,FBW,alpha):
#    t=arange(min(FIRlen,N-1)//2+1)*FBW+1e-10
#    Hf = FBW*(sin(pi*t*(1-alpha))+(4*alpha*t*cos(pi*t*(1+alpha))))/(pi*t*(1-(4*alpha*t)**2))
#    hf=concatenate([Hf,zeros(N-len(Hf)*2+1),flipud(Hf[1:])])
#    return(fft(hf/sum(hf)))

#far=lambda mu:array([1-mu,mu])   # Linear interpolation
far=lambda mu:array([mu*(-0.5+mu/2),-mu/2-mu**2/2+1,mu*(1.5-mu/2),-mu/2+mu**2/2]) #Farrow interpolation

def interp_farrow(xx,yy):
    y=[]
    yy=concatenate([[yy[0]],yy,[yy[-1]]*2])  # to extrapolate
    for x in xx:
        y.append(far(x%1)@yy[int(x):int(x)+4])
    return(array(y))



if True:
    fig1=figure(num=1,clear=True,facecolor='orange')
    x=0.07;y=0.9;xx=0.05;yy=0.04;
    Nfft_h         =TextBox(axes([x,y,xx,yy]),'FFT size 2^','11');y=y-yy
    Nifft_h        =TextBox(axes([x,y,xx,yy]),'IFFT size 2^','6');y=y-yy
    FIRlen_h       =TextBox(axes([x,y,xx,yy]),'FIRlen=','31');y=y-yy
    PolyWidth_h    =TextBox(axes([x,y,xx,yy]),'Poly Width=','4');y=y-yy
    OS_h           =TextBox(axes([x,y,xx,yy]),'Over Sample=','2');y=y-yy
    RollOff_h      =TextBox(axes([x,y,xx,yy]),'RollOff=','5%');y=y-yy
    RRCbeta_h      =TextBox(axes([x,y,xx,yy]),'RRCbeta=','0');y=y-yy
    #Nifft_h        =Slider(axes([0.04,0.35,0.02,0.6]),'IFFT',5,18,5,'%d',valstep=1,orientation='vertical')
    #OS_h           =Slider(axes([0.07,0.35,0.02,0.6]),'OS',1,OSfin,2,'%0.2f',valstep=0.1,orientation='vertical')
    Fc_h           =Slider(axes([0.13,0.35,0.02,0.6]),'Freq',-100,100,0,'%0.1f',valstep=0.2,orientation='vertical')
    #RollOff_h=RadioButtons(axes([0.01,0.21,.07,0.1]),('0%','5%','20%','35%'),active=1)
    OL_h     =RadioButtons(axes([0.08,0.21,.07,0.1]),('1/2','1/4','1/8','0'),active=1)
    Wtype_h  =RadioButtons(axes([0.01,0.05,.10,0.16]),('RECT OLdisc','BART OLadd','HANN OLadd','HANN OLdisc'),active=0)
    mode_h   =RadioButtons(axes([0.11,0.11,.05,0.1]),('RC-RRC','RC-FIR','POLY','PERFECT'),active=0)
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
        self.recalc(0)
    def read_GUI(self):
        self.Nfft=2**int(Nfft_h.text)
        self.Nifft=2**int(Nifft_h.text)
        self.OS=float(OS_h.text)
        self.OL=eval(OL_h.value_selected)
        self.alpha=int(RollOff_h.text[:-1])/100
        self.Wtype=Wtype_h.value_selected
        self.mode=mode_h.value_selected
        self.FIRlen=int(FIRlen_h.text)
        self.PolyWidth=int(PolyWidth_h.text)
        self.RRCbeta=float(RRCbeta_h.text)
    def standalone(self,Nfft,Nifft,OS,OL,alpha,Wtype,mode,FIRlen,PolyWidth,RRCbeta):
        
        self.Nfft=Nfft
        self.Nifft=Nifft
        self.OS=OS
        self.OL=OL
        self.alpha=alpha
        self.Wtype=Wtype
        self.mode=mode
        self.FIRlen=FIRlen
        self.PolyWidth=PolyWidth
        self.RRCbeta=RRCbeta
        return(self.IFFT_calc(False))
        
    def calculate(self,Fc,plotEn=True):
        """ Extracting variables 
        """
        Nbw=self.Nifft/self.OS*(1+self.alpha)
        OLwin=int(self.Nfft*self.OL)
        OLfft=int(self.Nfft*self.OL)
        OLifft=int(self.Nifft*self.OL)
        if 'OLdisc' in self.Wtype: OLwin//=2
        N=self.Nfft*os
        """ Creating Frequency envelope
        """
        n=int(round(Nbw*os/2*(1-self.alpha)/(1+self.alpha)))
        Fir=sqrt(rrc(0,(int(round(Nbw*os/2))-n),2*n))
        Fenv=concatenate([tile(Fir*10**(adj_ch_boost/20),(Nchan-1)//2),Fir,tile(Fir*10**(adj_ch_boost/20),(Nchan-1)//2),[0]])
        fe=arange(len(Fenv))/os-(len(Fenv)-1)/2/os
        """ Creating main signal (s)
        """
        t=arange(-N//2,N//2)
        s=exp(1j*2*pi*Fc/self.Nfft*t)*interp(Fc,fe,Fenv)
        """ Creating windows w: pre-FFT window wo: time envelope window 
        """
        if   'RECT' in self.Wtype:
            w=concatenate([zeros((N-self.Nfft)//2),ones(self.Nfft),zeros((N-self.Nfft)//2)])
            wo=concatenate([zeros((N-self.Nfft+OLfft)//2),ones(self.Nfft-OLfft),zeros((N-self.Nfft+OLfft)//2)])
        elif 'BART' in self.Wtype: 
            w=concatenate([zeros((N-self.Nfft)//2),
                           linspace(0,1,OLwin,endpoint=False),
                           ones(self.Nfft-2*OLwin),
                           linspace(1,0,OLwin,endpoint=False),
                           zeros((N-self.Nfft)//2)])
            wo=w
        elif 'HANN' in self.Wtype: 
            w=rrc((N-self.Nfft)//2,OLwin,self.Nfft-2*OLwin)
            wo=w
        """ Create overlapping FFT's s3
        """
        self.s3=zeros(span*self.Nifft-(span-1)*OLifft,'complex')
        f=arange(min(self.Nfft,self.Nifft*2)*os)-min(self.Nfft,self.Nifft*2)*os//2
        self.ww=zeros(N)
        for i in range(span):
            offset=(i-span//2)*(self.Nfft-OLfft)
            if plotEn:  L1d[i].set_data((t+offset)*self.Nifft/self.Nfft,wo[t+N//2])
            S=fftshift(fft(fftshift(w*roll(s,-offset))))/(self.Nfft-OLwin)
            if self.mode=='RC-RRC':  # select center BW only, fill with zeros outside
                n=int(ceil(Nbw/2)+trim)*2
                SS=S[arange(-n//2,n//2)*os+N//2]
            else:
                n=int(round((self.Nifft-Nbw)/2))
                SS=S[arange(-self.Nifft//2,self.Nifft//2)*os+N//2]*(rrc(0,n,self.Nifft-2*n))
            if i==span//2 and plotEn: # plot Frequency domain if center FFT
                L2a.set_data(f/os,real(S[f+N//2]))
                L2b.set_data(arange(-len(SS)//2,len(SS)//2),real(SS))
            ss=concatenate([zeros(self.Nifft//2-len(SS)//2),SS,zeros(self.Nifft//2-len(SS)//2)])
            ss=fftshift(ifft(fftshift(ss)))
            """ overlap and discard or add
            """
            if 'OLdisc' in self.Wtype:
                self.s3[arange(self.Nifft-OLifft)+i*(self.Nifft-OLifft)+OLifft//2]+=ss[arange(-(self.Nifft-OLifft)//2,(self.Nifft-OLifft)//2)+self.Nifft//2]*(self.Nifft-OLifft/2)
                self.ww[arange(-(self.Nfft-OLfft)//2,(self.Nfft-OLfft)//2)+offset+N//2]=1
            else:
                self.s3[arange(self.Nifft)+i*(self.Nifft-OLifft)]+=ss*(self.Nifft-OLifft)
                self.ww+=roll(wo,offset)
        """ Filter and interpolator
        """
        self.s3=fft(self.s3)
        if   self.mode=="RC-RRC": #filter first then oversample
             self.s4=self.RateConverter(self.s3)
             self.s4=self.s4*RRC1(len(self.s4),self.FIRlen,1/OSfin,self.alpha,self.RRCbeta)  #(1+self.alpha)
        if   self.mode=="RC-FIR": #filter first then oversample
             self.s4=self.RateConverter(self.s3)
             f2=interp(linspace(-1,1,len(self.s4)),linspace(-self.up,self.up,len(self.RCfir)),fftshift(self.up/(self.RCfir+1e-20)))
             self.s4=self.s4*RRC1(len(self.s4),self.FIRlen,1/OSfin,self.alpha,self.RRCbeta)*fftshift(f2)  #(1+self.alpha)
        elif self.mode=="POLY" or self.mode=="PERFECT": #oversample first then filter
             self.s4=self.RateConverter(self.s3)
        self.s4=ifft(self.s4)
        """ Creating reference signal (sr4)
        """
        sr=s*(abs(Fc)<Nbw/2)*self.ww*interp(Fc,fe,Fenv)
        i=(span*self.Nfft- (span-1)*OLfft)
        self.sr4=interp(linspace(-i/2,i/2,len(self.s4),endpoint=False),arange(-self.Nfft*os//2,self.Nfft*os//2),sr)
        """ Calculating MER
        """
        n=int(self.Nifft*(1-self.OL)/self.OS*OSfin/2*(span-2))
        r=arange(-n,n)+len(self.s4)//2
        err=abs(self.s4[r]-self.sr4[r])**2
        """ Plotting
        """
        t4=linspace(-len(self.s3)//2,len(self.s3)//2,len(self.s4))
        p4=mean(abs(self.sr4[r])**2)
        if plotEn:
            L1a.set_data(t4,real(self.sr4))
            L1b.set_data(t4,real(self.s4))
            f=arange(-len(self.s3)//2,len(self.s3)//2)
            L2e.set_data(f/len(self.s3)*self.Nifft,fftshift(self.RCfir)[f+len(self.RCfir)//2]/self.up)
            ax1.axis([-self.Nifft*os//2,self.Nifft*os//2,-1.5,1.5])
            textMER.set_text('MER spot=%0.1f'%(-10*log10(mean(err)+1e-20)))
            L2d.set_data(fe,Fenv)
            ax2.set_xticks([-self.Nifft//2,-Nbw/2,Nbw/2,self.Nifft//2])
            ax2.axis([-self.Nifft*0.6,self.Nifft*0.6,-0.5,1.5])
            fig1.canvas.draw()
        return(t4[r],err,p4)
    """_______________________________________________________________
    """

    def RateConverter(self,s):
        s=tile(s,self.up)*self.RCfir
        z=int(ceil(len(s)/self.dn))*self.dn-len(s)
        ss=concatenate([s[:len(s)//2],zeros(z),s[-len(s)//2:]])
        y1=sum(reshape(ss,[self.dn,-1]),0)/self.dn
        return(y1)
    def IFFT_calc(self,PlotEn=True):
        OLifft=int(self.Nifft*self.OL)
        (self.up,self.dn)=Fraction(OSfin/self.OS).limit_denominator().as_integer_ratio()
        n=int(span*self.Nifft-(span-1)*OLifft)*self.up
        if self.mode=="POLY":
            a=real(fftshift(ifft(RRC1(self.FIRlen,self.FIRlen,1/self.OS/self.PolyWidth,self.alpha,self.RRCbeta))))
            #b=interp(      arange(len(a)*self.up//self.PolyWidth)/self.up*self.PolyWidth,arange(len(a)),a)
            b=interp_farrow(arange(len(a)*self.up//self.PolyWidth)/self.up*self.PolyWidth,a)
            c=concatenate([zeros(n//2-len(b)//2),b,zeros(n//2-len(b)//2)])*self.PolyWidth
            self.RCfir=real(fft(fftshift(c)))
        elif self.mode=='PERFECT':
            self.RCfir=(RRC1(n,n,1/self.up/self.OS,self.alpha))*self.up
        else:
            fir=reshape(list(zip(*[far(i/self.up) for i in range(self.up,0,-1)])),[-1])
            self.RCfir=real(fft(fftshift(concatenate([zeros(n//2-len(fir)//2),fir,zeros(n//2-len(fir)//2)]))))  

        MERfreq=[]
        Fsweep=linspace(-self.Nifft,self.Nifft,NpointFsweep)
        mm=array([])
        pp=0
        for Fc in Fsweep:
            t,m,p=self.calculate(Fc,False)
            MERfreq.append(-10*log10(mean(m)+1e-20)/100)
            mm=mm+m if len(mm) else m
            pp+=p
            #print('%d %0.1f %0.2f'%(Fc,10*log10(mean(m)),pp))
        mer=-10*log10(mean(mm)/pp)
        if PlotEn:
            L1c.set_data(t,-10*log10(mm/pp+1e-20)/100)
            L2c.set_data(Fsweep,MERfreq)
            textMERavg.set_text('MER avg=%0.1fdB'%(mer))
        print('config: ',n,self.up,self.dn,mer,self.PolyWidth)
        return(mer)
    def recalc(self,a):
        self.read_GUI()
        self.IFFT_calc()
        self.calculate(Fc_h.val)
            
wf=waveform()

Fc_h.on_changed(wf.calculate)

#FIRlen_h.on_submit(wf.recalc)
#PolyWidth_h.on_submit(wf.recalc)
#Nfft_h.on_submit(wf.recalc)
#Nifft_h.on_submit(wf.recalc)
#OS_h.on_submit(wf.recalc)
#
#RollOff_h.on_submit(wf.recalc)
#OL_h.on_clicked(wf.recalc)
Wtype_h.on_clicked(wf.recalc)
#mode_h.on_clicked(wf.recalc)

show()


