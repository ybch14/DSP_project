close all;clear;clc;
f1=324;f2=310;f3=324.4;
Fs=1000;
t=0:1/Fs:6;NFFT=1024*2^nextpow2(length(t));
f=Fs/2*linspace(0,1,NFFT/2+1);
f_e=zeros(100,3);
for time=1:100
    y=5*exp(1i*2*pi*f1*t)+5*exp(1i*2*pi*f2*t)+5*exp(1i*2*pi*f3*t)+10*randn(1,length(t));
    F=fft(y,NFFT)/length(t);
    F=2*abs(F(1:NFFT/2+1));
    % plot(f,F);
    [~,f_e(time,:)]=findpeaks(F,f,'MinPeakHeight',max(F)*0.5);
end
E=mean(f_e)
MSE=mean((f_e-repmat(mean(f_e),100,1)).^2)