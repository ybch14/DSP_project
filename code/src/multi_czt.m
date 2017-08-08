close all;clear;clc;
Fs=1000;
f1=324;f2=310;f3=324.12;
t=0:1/Fs:6;
NFFT=1024*2^nextpow2(length(t));
f=Fs/2*linspace(0,1,NFFT/2+1);
y=5*exp(1i*2*pi*f1*t)+5*exp(1i*2*pi*f2*t)+5*exp(1i*2*pi*f3*t)+2.5*randn(1,length(t));
F=fft(y,NFFT)/length(t);
F=2*abs(F(1:NFFT/2+1));
[~,f_pre]=findpeaks(F,f,'MinPeakHeight',max(F)*0.5);sort(f_pre);
f_low=f_pre(1)*0.98;f_high=f_pre(3)*1.02;


delta_omega=20000;

omega=exp(-1i*2*pi*(f_high-f_low)/(Fs*delta_omega));
omega_start=exp(1i*2*pi*f_low/Fs);
h=0:1:delta_omega-1;
f=(f_high-f_low)/delta_omega*h+f_low;
f_e=zeros(100,3);
for time=1:100
    y=5*exp(1i*2*pi*f1*t)+5*exp(1i*2*pi*f2*t)+5*exp(1i*2*pi*f3*t)+10*randn(1,length(t));
    x=czt(y,delta_omega,omega,omega_start);
    [~,locs]=findpeaks(abs(x),f,'MinPeakHeight',max(abs(x))/2);
    sort(locs);f_e(time,:)=locs;
end
figure;
plot(f,abs(x));
figure;
plot(1:time,f_e(:,1)-f2,'r');
hold on;
plot(1:time,f_e(:,2)-f1,'g');
hold on;
plot(1:time,f_e(:,3)-f3,'b');
E=mean(f_e)
MSE=mean((f_e-repmat([f2,f1,f3],100,1)).^2)


