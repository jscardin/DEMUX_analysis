from pylab import *
from matplotlib import animation
from matplotlib.widgets import Slider, Button, RadioButtons
from quicktions import Fraction

os=8   # oversampling for the entire simulation
OSfin=4  #final oversampling
span=7 # number of FFT frame, must be odd number
FIRlen=512 # FIR length

def rrc(N0,N01,N1):
    return(concatenate([zeros(N0),
                        cos(linspace(-pi,0,N01,endpoint=False))/2+0.5,
                        ones(N1),
                        cos(linspace(0,pi,N01,endpoint=False))/2+0.5,
                        zeros(N0)]))

def RRC1(N,FBW,alpha):
    ''' returns zero phase RRC impulse response'''
    n=int(round(N/2*FBW))
    z=int(round(n*alpha))
    H=concatenate([ones(n-z),sqrt((cos(linspace(0,pi,2*z,endpoint=False))+1)/2)])
    H=concatenate([H,zeros(N-2*len(H)+1),flipud(H[1:])])
    return(H)



#fig2,ax=subplots(1,num=2,clear=True)

#fig,(a,ax1,ax2)=subplots(1,3,num=1,clear=True)
fig=figure(num=1,clear=True,facecolor='orange')

Nfft_h         =Slider(axes([0.01,0.35,0.02,0.6]),'FFT',5,18,10,'%d',valstep=1,orientation='vertical')
Nifft_h        =Slider(axes([0.04,0.35,0.02,0.6]),'IFFT',5,18,7,'%d',valstep=1,orientation='vertical')
OS_h           =Slider(axes([0.07,0.35,0.02,0.6]),'OS',1,OSfin,2,'%0.2f',valstep=0.25,orientation='vertical')
Fc_h           =Slider(axes([0.10,0.35,0.02,0.6]),'Freq',-100,100,0,'%0.1f',valstep=0.2,orientation='vertical')
Toff_h         =Slider(axes([0.13,0.35,0.02,0.6]),'Time',-200,200,0,'%d',orientation='vertical',valstep=1)
RollOff_h=RadioButtons(axes([0.01,0.21,.07,0.1]),('0%','5%','20%','35%'),active=1)
OL_h     =RadioButtons(axes([0.08,0.21,.07,0.1]),('1/2','1/4','1/8','0'),active=1)
Wtype_h  =RadioButtons(axes([0.01,0.01,.10,0.2]),('RECT OLdisc','BART OLadd','HANN OLadd','HANN OLdisc'),active=2)
order_h   =RadioButtons(axes([0.11,0.11,.05,0.1]),('RC-FIR','FIR-RC','MASK-RC'),active=1)
RCtype_h   =RadioButtons(axes([0.11,0.01,.05,0.1]),('Linear','Farrow','FLAT'),active=2)


ax1=subplot2grid((8,10),(0,1),rowspan=4,colspan=9,title='Time Domain Samples')
L1a,=ax1.plot([],[],'y.-',lw=1,ms=2,label='real(Input Signal)')
L1b,=ax1.plot([],[],'c.-' ,lw=1,ms=2,label='real(DEMUX output)')
L1c,=ax1.plot([],[],'b',lw=2,label='EVM')
L1d=[ax1.plot([],[],'r--',lw=1,label='FFT blocks')[0]]

for i in range(span-1):
    L1d.append(ax1.plot([],[],'r--',lw=1)[0])
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_xlabel('IFFT TIME SAMPLES')


ax2=subplot2grid((8,10),(5,1),rowspan=3,colspan=9,title='Freq Domain Tones')
L2a,=ax2.plot([],[],'y',lw=4,label='real(Input Signal)')
L2b,=ax2.plot([],[],'c.',ms=3,label='real(IFFT input)')
L2c,=ax2.plot([],[],'r',label='Channel Envelope')
L2d,=ax2.plot([],[],label='Interpolator Freq resp')
L2e,=ax2.plot([],[],'b',lw=2,label='EVM (%)')
ax2.legend()
ax2.grid()
ax2.set_xlabel('IFFT FREQ TONES')



textMER=ax1.text(0,-1.5,'asdf',size='30',ha='center')


class waveform():
    def __init__(self):
        pass

 

    def calculate(self,Fc,plotEn=True):
        #______ Extracting variables
        #Fcs=Fc_h.val
        Nfft=int(2**Nfft_h.val)
        Nifft=int(2**Nifft_h.val)
        Nbw=Nifft/OS_h.val
        OLwin=int(Nfft*eval(OL_h.value_selected))
        OLfft=int(Nfft*eval(OL_h.value_selected))
        OLifft=int(Nifft*eval(OL_h.value_selected))
        if 'OLdisc' in Wtype_h.value_selected: OLwin//=2
        ro=[a[0] for a in ((0,"0%"),(.05,"5%"),(.2,"20%"),(.35,"35%")) if a[1]==RollOff_h.value_selected][0]
        N=Nfft*os
        #______ Creating Frequency envelope
        n=int(round(Nbw*os/2*(1-ro)/(1+ro)))
        Fir=sqrt(rrc(0,(int(round(Nbw*os/2))-n),2*n))
        Fenv=concatenate([tile(Fir,span),[0]])
        fe=arange(len(Fenv))/os-(len(Fenv)-1)/2/os
        #______ Creating main signal (s)
        t=arange(-N//2,N//2)
        s=exp(1j*2*pi*Fc/Nfft*(t+Toff_h.val))*interp(Fc,fe,Fenv)
        #______ Creating windows w: pre-FFT window wo: time envelope window 
        if 'RECT' in Wtype_h.value_selected:
            w=concatenate([zeros((N-Nfft)//2),ones(Nfft),zeros((N-Nfft)//2)])
            wo=concatenate([zeros((N-Nfft+OLfft)//2),ones(Nfft-OLfft),zeros((N-Nfft+OLfft)//2)])
        elif 'BART' in Wtype_h.value_selected: 
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
            if plotEn:  L1d[i].set_data((t+offset)*Nifft/Nfft,wo[t+N//2])
            S=fftshift(fft(fftshift(w*roll(s,-offset))))/(Nfft-OLwin)
            if order_h.value_selected=="RC-FIR":  # select center BW only, fill with zeros outside
                n=int(round(Nbw/2)*2)
                SS=S[arange(-n//2,n//2)*os+N//2]
            elif order_h.value_selected=="FIR-RC":#Select whole range
                n=int(round((Nifft-Nbw)/2))
                SS=S[arange(-Nifft//2,Nifft//2)*os+N//2]*(rrc(0,n,Nifft-2*n))
            elif order_h.value_selected=="MASK-RC":#apply RRC filter and RC compensation
                n=int(round(Nbw/2*(1-ro)/(1+ro)))
                nn=int(round(Nbw/2)*2)
                f1=sqrt(rrc(0,(int(round(Nbw/2))-n),2*n))
                f2=interp(arange(-nn//2,nn//2),linspace(-Nifft//2*self.up,Nifft//2*self.up,len(self.RCfir)),fftshift(self.up/self.RCfir))
                SS=S[arange(-nn//2,nn//2)*os+N//2]*f1*f2
            if i==span//2 and plotEn: # plot Frequency domain if center FFT
                L2a.set_data(f/os,real(S[f+N//2]))
                L2b.set_data(arange(-len(SS)//2,len(SS)//2),real(SS))
            ss=concatenate([zeros(Nifft//2-len(SS)//2),SS,zeros(Nifft//2-len(SS)//2)])
            ss=fftshift(ifft(fftshift(ss)))
            #______ overlap and discard or add
            if 'OLdisc' in Wtype_h.value_selected:
                self.s3[arange(Nifft-OLifft)+i*(Nifft-OLifft)+OLifft//2]+=ss[arange(-(Nifft-OLifft)//2,(Nifft-OLifft)//2)+Nifft//2]*(Nifft-OLifft/2)
                self.ww[arange(-(Nfft-OLfft)//2,(Nfft-OLfft)//2)+offset+N//2]=1
            else:
                self.s3[arange(Nifft)+i*(Nifft-OLifft)]+=ss*(Nifft-OLifft)
                self.ww+=roll(wo,offset)
        #______ Filter and interpolator
        self.s3=fft(self.s3)
        if order_h.value_selected=="RC-FIR": #filter first then oversample
            self.s4=self.RateConverter(self.s3)
            self.s4=self.s4*RRC1(len(self.s4),1/OSfin/(1+ro),ro)
        elif order_h.value_selected=="FIR-RC": #oversample first then filter
            self.s3=self.s3*RRC1(len(self.s3),Nbw/Nifft/(1+ro),ro)
            self.s4=self.RateConverter(self.s3)
        else: #oversample first then filter
            self.s4=self.RateConverter(self.s3)
        self.s4=ifft(self.s4)
        #______ Creating reference signal (sr4)
        sr=s*(abs(Fc)<Nbw/2)*self.ww*interp(Fc,fe,Fenv)
        i=(span*Nfft- (span-1)*OLfft)
        self.sr4=interp(linspace(-i/2,i/2,len(self.s4),endpoint=False),arange(-Nfft*os//2,Nfft*os//2),sr)
        #______ Calculating MER
        n=int(Nifft*OSfin/(Nifft/Nbw)*1)
        r=arange(-n,n)+len(self.s4)//2
        err=abs(self.s4[r]-self.sr4[r])**2
        #______ Plotting
        t4=arange(0,len(self.s3),(Nifft/Nbw)/OSfin)-len(self.s3)//2
        p4=mean(abs(self.sr4[r])**2)
        if plotEn:
            L1a.set_data(t4,real(self.sr4))
            L1b.set_data(t4,real(self.s4))
            f=arange(-len(self.s3)//2,len(self.s3)//2)
            L2d.set_data(f/len(self.s3)*Nifft,fftshift(self.RCfir)[f+len(self.RCfir)//2]/self.up)
            ax1.axis([-Nifft*os//2,Nifft*os//2,-1.5,1.5])
            textMER.set_text('MER avg=%0.1fdB   MER spot=%0.1f'%(wf.MERtotal,-10*log10(mean(err)+1e-20)))
            L2c.set_data(fe,Fenv)
            ax2.set_xticks([-Nifft//2,-Nbw/2,Nbw/2,Nifft//2])
            ax2.axis([-Nifft*0.6,Nifft*0.6,-0.5,1.5])
            fig.canvas.draw()
            #fig.canvas.flush_events()
        return(t4[r],err,p4)
    #_____________________________________________________________________
    def calc(self,a):
        self.calculate(Fc_h.val)

    def IFFT_calc(self,a):
        Nifft=int(2**Nifft_h.val)
        OLifft=int(Nifft*eval(OL_h.value_selected))
        n=span*Nifft-(span-1)*OLifft
        (self.up,self.dn)=Fraction(OSfin/OS_h.val).limit_denominator().as_integer_ratio()
        if RCtype_h.value_selected=="FLAT":
            self.RCfir=fftshift(concatenate([zeros(n*self.up//2-n//2),ones(n),zeros(n*self.up//2-n//2)]))*self.up     
        else:
            if RCtype_h.value_selected=="Linear":
                far=lambda mu:array([1-mu,mu])
            else: #RCtype_h.value_selected=="Farrow":
                far=lambda mu:array([mu*(-0.5+mu/2),-mu/2-mu**2/2+1,mu*(1.5-mu/2),-mu/2+mu**2/2])
            fir=reshape(list(zip(*[far(i/self.up) for i in range(self.up,0,-1)])),[-1])
            self.RCfir=real(fft(fftshift(concatenate([zeros(n*self.up//2-len(fir)//2),fir,zeros(n*self.up//2-len(fir)//2)]))))  
        self.SweepFreq()
        self.calculate(Fc_h.val)
    def RateConverter(self,s):
        s=tile(s,self.up)*self.RCfir
        z=int(ceil(len(s)/self.dn))*self.dn-len(s)
        ss=concatenate([s[:len(s)//2],zeros(z),s[-len(s)//2:]])
        y1=sum(reshape(ss,[self.dn,-1]),0)/self.dn
        return(y1)
    def SweepFreq(self):
        MERfreq=[]
        N=500
        F=linspace(-2**Nifft_h.val,2**Nifft_h.val,N)
        mm=array([])
        pp=0
        for f in F:
            t,m,p=self.calculate(f,False)
            MERfreq.append(-10*log10(mean(m)+1e-20)/100)
            mm=mm+m if len(mm) else m
            pp+=p
            #print('%d %0.1f %0.2f'%(f,10*log10(mean(m)),pp))
        L1c.set_data(t,-10*log10(mm/pp+1e-20)/100)
        self.MERtotal=-10*log10(mean(mm)/pp)
        L2e.set_data(F,MERfreq)
    def SweepOS(self):
        OSrange=linspace(1,4,13)
        figure(2);clf()
        for ni in range(5,9):
            print('____ Nifft=',2**ni)
            Nifft_h.set_val(ni)
            mer=[]
            for o in OSrange:
                OS_h.set_val(o)
                self.SweepFreq()
                mer.append(self.MERtotal)
                print(o,self.MERtotal)
            plot(OSrange,mer)
        grid()
        legend()
            
wf=waveform()


wf.IFFT_calc(0)

Fc_h.on_changed(wf.calculate)

Toff_h.on_changed(wf.calc)
Nfft_h.on_changed(wf.IFFT_calc)

Nifft_h.on_changed(wf.IFFT_calc)
OS_h.on_changed(wf.IFFT_calc)
RollOff_h.on_clicked(wf.IFFT_calc)
OL_h.on_clicked(wf.IFFT_calc)
Wtype_h.on_clicked(wf.IFFT_calc)
order_h.on_clicked(wf.IFFT_calc)
RCtype_h.on_clicked(wf.IFFT_calc)

show()