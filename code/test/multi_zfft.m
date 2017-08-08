close all;clear;clc;
f1=324;f2=310;f3=324.3;
f_low=300;f_high=330;
Fs=1000;
t=0:1/Fs:6;
f_e=zeros(100,3);
for time=1:100
    y=5*exp(1i*2*pi*f1*t)+5*exp(1i*2*pi*f2*t)+5*exp(1i*2*pi*f3*t)+10*randn(1,length(t));
    [f,F]=my_zfft(y,f_low,f_high,Fs);
    [~,locs]=findpeaks(F,f,'MinPeakHeight',max(F)/2);
    sort(locs);f_e(time,:)=locs;
end
plot(f,F);
E=mean(f_e)
MSE=mean((f_e-repmat([f2,f1,f3],100,1)).^2)