from pylab import *
from matplotlib import animation
from matplotlib.widgets import Slider, Button, RadioButtons


Fs2=50
Nframe=64
overlap=1/4
OL=0.1
Fc=-0
anim_enable=False


#fig,(a,ax1,ax2)=subplots(1,3,num=1,clear=True)
fig=figure(num=1,clear=True)

Nfft_h=Slider(axes([0.01,0.3,0.02,0.6]),'Nfft',5,20,16,'%d',valstep=1,orientation='vertical')
Nifft_h=Slider(axes([0.04,0.3,0.02,0.6]),'Nifft',5,20,8,'%d',valstep=1,orientation='vertical')
Nidft_h=Slider(axes([0.07,0.3,0.02,0.6]),'Nidft',1,100,50,'%d',orientation='vertical',valstep=1)
Fc_h=Slider(axes([0.1,0.3,0.02,0.6]),'Fc',-32,32,0,'%0.1f',orientation='vertical')
Toff_h=Slider(axes([0.13,0.3,0.02,0.6]),'Toff',-200,200,0,'%d',orientation='vertical',valstep=1)
Stim_h=RadioButtons(axes([0.11,0.01,.05,0.1]),('CW','Impulse'),active=0)
OL_h=RadioButtons(axes([0.01,0.1,.1,0.1]),('1/2','1/4','1/8','0'),active=1)
Wtype_h=RadioButtons(axes([0.01,0.01,.1,0.1]),('RECT','BARTLETT','HANN'),active=2)


ax1=subplot2grid((8,10),(0,1),rowspan=4,colspan=9,title='Time Domain \n(Real part only, imag not shown)')
L1InReal,=ax1.plot([],[],'y',lw=5,label='real(DEMUX input)')
#L1OutReal.append(ax1.plot([],[],'r.-',lw=1,label='DEMUX output (Real)')[0])
#L1OutReal.append(ax1.plot([],[],'g.-',lw=1,label='DEMUX output (Real)')[0])
#L1OutReal.append(ax1.plot([],[],'r.-',lw=1,label='DEMUX output (Real)')[0])
#L1imag,=ax1.plot([],[],lw=1,label='Input Waveform')
L1OutReal,=ax1.plot([],[],'b',lw=1,label='real(DEMUX output)')
L1win1,=ax1.plot([],[],'r',lw=1,label='FFT block 1')
L1win2,=ax1.plot([],[],'r',lw=1,label='FFT block 2')
L1win3,=ax1.plot([],[],'r',lw=1,label='FFT block 3')
ax1.legend(loc='upper right')



ax2=subplot2grid((8,10),(5,1),rowspan=3,colspan=9)
L2fftout,=ax2.plot([],[],lw=1,label='mag(FFT out)')
L2fftoutreal,=ax2.plot([],[],lw=1,label='real(FFT out)')
L2fftpoints,=ax2.plot([],[],'.',label='real(IFFT input)')
ax2.legend()

#ax3=subplot2grid((8,10),(0,1),rowspan=1,colspan=3)
#ax3.imshow(img)
#ax3.axis('off')

#ax1.set_title('Input Signal')
ax1.grid()
ax1.set_xlabel('Time (samples)')


ax2.set_title('Input Signal Spectrum (infinite average)')
ax2.grid()
ax2.set_xlabel('FREQ')



textMER=ax1.text(0,-1.5,'asdf')


class waveform():
    def __init__(self):
        self.calculate(1)


    def calculate(self,a):
        Nfft=int(2**Nfft_h.val)
        Nifft=int(2**Nifft_h.val)
        Nidft=int(round(Nidft_h.val/100*Nifft))
        OL=int(Nfft*eval(OL_h.value_selected))
        ov=[OL-Nfft,0,Nfft-OL]  #index of the center of the 3 FFT'self.s
        N=int(Nframe*Nifft*(1-overlap))
        N=2**20
        print(N)
        t=arange(-N//2,N//2)
        if Stim_h.value_selected=='CW':
            self.s=exp(1j*2*pi*Fc_h.val/Nfft*(t+Toff_h.val))
            ss=self.s*(abs(Fc_h.val)<Nidft//2)
            
        else:
            self.s=zeros(N)
            self.s[int(Toff_h.val)+N//2]=Nifft/Nidft
            ss=rfft(self.s)
            ss[round(Nidft/Nifft*N/2):]=0
            ss=irfft(ss)

        if Wtype_h.value_selected=='RECT':
            w=zeros(N)
            w[arange(-(Nfft-OL)//2,(Nfft-OL)//2)+N//2]=1
        elif Wtype_h.value_selected=='BARTLETT': 
            w=concatenate([zeros((N-Nfft)//2),linspace(0,1,OL,endpoint=False),ones(Nfft-2*OL),linspace(1,0,OL,endpoint=False),zeros((N-Nfft)//2)])
        elif Wtype_h.value_selected=='HANN': 
            w=concatenate([zeros((N-Nfft)//2),cos(linspace(-pi,0,OL,endpoint=False))/2+0.5,ones(Nfft-2*OL),cos(linspace(0,pi,OL,endpoint=False))/2+0.5,zeros((N-Nfft)//2)])

        ww=zeros(N)
        for i in range(3):
            ww+=roll(w,ov[i])
        L1InReal.set_data(t,self.s*ww)
        #L1imag.set_data(t,imag(self.s))
        L1win1.set_data(t-ov[0],w[t+N//2])
        L1win2.set_data(t,w[t+N//2])
        L1win3.set_data(t+ov[0],w[t+N//2])

#        s3=zeros(3*Nfft-2*OL,'complex')
#        f=arange(-Nifft//2,Nifft//2)
#        for i in range(3):
#            if Wtype_h.value_selected=='RECT':
#                www=concatenate([zeros((N-Nfft)//2),ones(Nfft),zeros((N-Nfft)//2)])
#                S=fftshift(fft(fftshift(www*roll(self.s,-ov[i]))))
#            else:
#                S=fftshift(fft(fftshift(w*roll(self.s,-ov[i]))))
#            #S=roll(S,int(Fc_h.val*N/Nifft))
#            #S=convolve(S,ones(int(OL*N)+1),'same')#/(int(OL*N)+1)
#            if i==1: 
#                L2fftout.set_data(    f, abs(S[f+N//2])/Nifft)
#                L2fftoutreal.set_data(f,real(S[f+N//2])/Nifft)
#            S[0:(Nifft-Nidft)//2]=0
#            S[(Nifft+Nidft)//2:]=0
#            if i==1: L2fftpoints.set_data(arange(-Nidft//2,Nidft//2)+1,real(S[arange(-Nidft//2,Nidft//2)+1+Nifft//2])/Nifft)
#            print(len(S),Nifft,Nidft,len(w),Toff_h.val)
#            S=fftshift(ifft(fftshift(S)))
#            #L1OutReal[i+1].set_data(arange(-Nifft//2,Nifft//2)+ov[i],real(SS))
#            if Wtype_h.value_selected=='RECT':                
#                s3[arange(Nfft-OL)+i*(Nfft-OL)+OL//2]+=S[OL//2:Nfft-OL//2]
#            else:
#                s3[arange(Nfft)+i*(Nfft-OL)]+=S[arange(-Nfft//2,Nfft//2)+ov[i]+N//2]
#        L1OutReal.set_data(arange(OL-(Nifft*3)//2,(Nifft*3)//2-OL),s3)
#        i=arange(-len(s3)//2,len(s3)//2)+N//2
#        #textMER.set_text('MER est.=%0.1f'%(10*log10(sum(abs(s3-ss[i]*ww[i])**2)/sum(abs(ww[i])**2))))
#        
#        ax1.axis([-Nifft*2,Nifft*2,-1.5,1.5])
#        a=array([-Nifft*3//2+OL,-Nifft*3//2+2*OL,-Nifft*3//2+OL,-Nifft//2,-Nifft//2+OL])
#        ax1.set_xticks(concatenate([a,-a]))
#        ax2.set_xticks([-Nifft//2,-Nidft//2+1,Nidft//2,Nifft//2])
#        ax2.axis([-Nifft/2,Nifft/2,-1,2])
        fig.canvas.draw()
        fig.canvas.flush_events()
        

wf=waveform()

fig.canvas.toolbar.update()

Nifft_h.on_changed(wf.calculate)
Nidft_h.on_changed(wf.calculate)
Fc_h.on_changed(wf.calculate)
Toff_h.on_changed(wf.calculate)
OL_h.on_clicked(wf.calculate)
Wtype_h.on_clicked(wf.calculate)

