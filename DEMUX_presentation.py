from pylab import *
from matplotlib import animation
from matplotlib.widgets import Slider, Button, RadioButtons

os=4


#fig,(a,ax1,ax2)=subplots(1,3,num=1,clear=True)
fig=figure(num=1,clear=True)

Nfft_h=Slider(axes([0.01,0.3,0.02,0.6]),'FFT',5,18,12,'%d',valstep=1,orientation='vertical')
Nifft_h=Slider(axes([0.03,0.3,0.02,0.6]),'IFFT',5,18,7,'%d',valstep=1,orientation='vertical')
Nidft_h=Slider(axes([0.05,0.3,0.02,0.6]),'IDFT',0,1,0.5,'%0.2f',orientation='vertical')
Nfir_h=Slider(axes([0.07,0.3,0.02,0.6]),'FIR',0,1,1,'%0.2f',orientation='vertical')
Fc_h=Slider(axes([0.1,0.3,0.02,0.6]),'Freq',-100,100,0,'%0.1f',valstep=0.2,orientation='vertical')
Toff_h=Slider(axes([0.12,0.3,0.02,0.6]),'Time',-200,200,0,'%d',orientation='vertical',valstep=1)
RollOff_h=RadioButtons(axes([0.08,0.01,.07,0.1]),('0%','5%','20%','35%'),active=3)
OLwin_h=RadioButtons(axes([0.01,0.11,.07,0.1]),('1/2','1/4','1/8','0'),active=1)
OL_h=RadioButtons(axes([0.08,0.11,.07,0.1]),('1/2','1/4','1/8','0'),active=1)
Wtype_h=RadioButtons(axes([0.01,0.01,.07,0.1]),('RECT','BARTLETT','HANN'),active=0)
FIRlen=512


ax1=subplot2grid((8,10),(0,1),rowspan=4,colspan=9,title='Time Domain Samples')
L1InReal,=ax1.plot([],[],'y',lw=4,label='real(Input Signal)')
L1OutReal,=ax1.plot([],[],'b.',lw=1,label='real(DEMUX output)')
L1win=[]
L1win.append(ax1.plot([],[],'r',lw=1,label='FFT blocks')[0])
L1win.append(ax1.plot([],[],'r',lw=1)[0])
L1win.append(ax1.plot([],[],'r',lw=1)[0])
ax1.legend(loc='upper right')
ax1.grid()
ax1.set_xlabel('TIME SAMPLES')

def rrc(N0,N01,N1):
    return(concatenate([zeros(N0),
                        cos(linspace(-pi,0,N01,endpoint=False))/2+0.5,
                        ones(N1),
                        cos(linspace(0,pi,N01,endpoint=False))/2+0.5,
                        zeros(N0)]))



ax2=subplot2grid((8,10),(5,1),rowspan=3,colspan=9,title='Freq Domain Tones')
#L2fftout,=ax2.plot([],[],lw=1,label='mag(FFT out)')
L2fftoutreal,=ax2.plot([],[],'y',lw=4,label='real(Input Signal)')
L2fftpoints,=ax2.plot([],[],'b.',label='real(IFFT input)')
L2env,=ax2.plot([],[],'r',label='Channel Envelope')
ax2.legend()
ax2.grid()
ax2.set_xlabel('FREQ TONES')



textMER=ax1.text(0,-1.5,'asdf')


class waveform():
    def __init__(self):
        self.calculate(1)


    def calculate(self,a):
        Nfft=int(2**Nfft_h.val)
        Nifft=int(2**Nifft_h.val)
        Nidft=int(round(Nidft_h.val*Nifft))
        Nfir=Nfir_h.val*Nifft
        OLwin=int(Nfft*eval(OL_h.value_selected))
        OLfft=int(Nfft*eval(OL_h.value_selected))
        OLifft=int(Nifft*eval(OL_h.value_selected))
        Fc=Fc_h.val*Nifft/160
        ro=[a[0] for a in ((0,"0%"),(.05,"5%"),(.2,"20%"),(.35,"35%")) if a[1]==RollOff_h.value_selected][0]
        ov=[OLfft-Nfft,0,Nfft-OLfft]  #index of the center of the 3 FFT's
        N=Nfft*os
        n=os*min(Nfir,Nidft)
        Fir=rrc(0,int(round(n*ro/2)),int(round(n*(1-ro))))
        Fenv=concatenate([Fir,Fir,Fir])
        f=arange(len(Fenv))/os-len(Fenv)/2/os
        L2env.set_data(f,Fenv)
        t=arange(-N//2,N//2)
        s=exp(1j*2*pi*Fc/Nfft*(t+Toff_h.val))*interp(Fc,f,Fenv)
        sr=s*(abs(Fc)<min(Nidft,Nfir)/2)
        if Wtype_h.value_selected=='RECT':
            w=zeros(N)
            w[arange(-(Nfft-OLfft)//2,(Nfft-OLfft)//2)+N//2]=1
        elif Wtype_h.value_selected=='BARTLETT': 
            w=concatenate([zeros((N-Nfft)//2),
                           linspace(0,1,OLwin,endpoint=False),
                           ones(Nfft-2*OLwin),
                           linspace(1,0,OLwin,endpoint=False),
                           zeros((N-Nfft)//2)])
        elif Wtype_h.value_selected=='HANN': 
            w=rrc((N-Nfft)//2,OLwin,Nfft-2*OLwin)        
        ww=zeros(N)
        for i in range(3):
            ww+=roll(w,ov[i])
            L1win[i].set_data(t+ov[i],w[t+N//2])
        L1InReal.set_data(t,real(sr)*ww)
        
        
        s3=zeros(3*Nifft-2*OLifft,'complex')
        f=arange(Nifft*os*2)-Nifft*os
        ff=arange(Nidft)-Nidft//2
        #L2fir.set_data(f,real(fftshift(FIR4)))
        for i in range(3):
            if Wtype_h.value_selected=='RECT':
                www=concatenate([zeros((N-Nfft)//2),ones(Nfft),zeros((N-Nfft)//2)])
                S=fftshift(fft(fftshift(www*roll(s,-ov[i]))))/Nfft
            else:
                S=fftshift(fft(fftshift(w*roll(s,-ov[i]))))/(Nfft-OLfft)


            if i==1: 
        #                L2fftout.set_data(    f/os, abs(S[f+N//2]))
                L2fftoutreal.set_data(f/os,real(S[f+N//2]))
            SS=zeros(Nifft,'complex')
            SS[ff+Nifft//2]=S[ff*os+N//2]
            #S*=FIR4
            if i==1: L2fftpoints.set_data(ff,real(SS[ff+Nifft//2]))
            ss=fftshift(ifft(fftshift(SS)))
            if Wtype_h.value_selected=='RECT' or OLwin<OLfft:                
                s3[arange(Nifft-OLifft)+i*(Nifft-OLifft)+OLifft//2]+=ss[arange(-(Nifft-OLifft)//2,(Nifft-OLifft)//2)+Nifft//2]*Nifft
            else:
                s3[arange(Nifft)+i*(Nifft-OLifft)]+=ss*(Nifft-OLifft)
        if Nfir<Nidft:
            a=rrc(round(FIRlen*(3-ro)/2/os),round(ro*FIRlen/os),round((1-ro)*FIRlen//os))
            a=concatenate([zeros(os*Nifft//2-len(Fir)//2),Fir,zeros(os*Nifft//2-len(Fir)//2)])
            print(len(s3),len(a))
            s3=convolve(s3,fftshift(ifft(fftshift(a))),'same')
            print(len(s3))
        L1OutReal.set_data((arange(len(s3))-len(s3)//2)*Nfft//Nifft,real(s3))
        i=arange(-len(s3)//2,len(s3)//2)*Nfft//Nifft+N//2
        textMER.set_text('MER est.=%0.1f'%(10*log10(sum(abs(s3-sr[i]*ww[i])**2)/sum(abs(ww[i])**2)+1e-20)))
        ax1.axis([-Nfft*2,Nfft*2,-1.5,1.5])
        a=array([-Nfft*3//2+OLifft,-Nfft*3//2+2*OLifft,-Nfft*3//2+OLifft,-Nfft//2,-Nfft//2+OLifft])
        #ax1.set_xticks(concatenate([a,-a]))
        ax2.set_xticks([-Nifft//2,-Nidft//2,Nidft//2,Nifft//2])
        ax2.axis([-Nifft*0.6,Nifft*0.6,-0.5,1.5])
        fig.canvas.draw()
        #fig.canvas.flush_events()
        

wf=waveform()

fig.canvas.toolbar.update()

Nfft_h.on_changed(wf.calculate)
Nifft_h.on_changed(wf.calculate)
Nidft_h.on_changed(wf.calculate)
Fc_h.on_changed(wf.calculate)
Toff_h.on_changed(wf.calculate)
Nfir_h.on_changed(wf.calculate)
RollOff_h.on_clicked(wf.calculate)
OL_h.on_clicked(wf.calculate)
Wtype_h.on_clicked(wf.calculate)

