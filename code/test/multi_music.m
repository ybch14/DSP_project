close all;clear;clc;
f1=324;f2=310;f3=326;
Fs=1000;
t=0:1/Fs:6;
f_e=zeros(100,3);
for time=1:100
    y=5*exp(1i*2*pi*f1*t)+5*exp(1i*2*pi*f2*t)+5*exp(1i*2*pi*f3*t)+2.5*randn(1,length(t));
    [omega,w,P]=my_music(y,3,100,5000);
    f_e(time,:)=sort(omega*Fs/(2*pi));
end
plot(w*Fs/(2*pi),P);
E=mean(f_e)
MSE=mean((f_e-repmat([f2,f1,f3],100,1)).^2)
