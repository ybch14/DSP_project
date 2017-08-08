function [f,F]=my_zfft(y,fi,fa,fs)
% y is signal
% fi is the start frequency
% fa is the end frequency
% fs is the sample rate of the signal
% w is the digital frequency
% F is the output spectrum
fe=(fi+fa)/2;
N=length(y);
r=0:N-1;
b=2*pi*fe.*r./fs;
temp1=y.*exp(-1i.*b);
Bw=fa-fi;
B=fir1(32,Bw/fs);
temp2=filter(B,1,temp1);
c=temp2(1:floor(fs/Bw/2):N);
Nc=length(c);
f=linspace(fi,fa,Nc);
F=abs(fft(c))./Nc*2;
F=circshift(F,[0,floor(Nc/2)]);
end