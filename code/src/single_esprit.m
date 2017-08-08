close all;clear;clc;
f=324.2333;
Fs=1000;
t=0:1/Fs:6;
f_e=zeros(1,100);
for time=1:100
    y=5*exp(1i*2*pi*f*t)+5*randn(1,length(t));
    omega=my_esprit(y,1,100);
    f_e(time)=omega*Fs/(2*pi);
end
figure;
plot(1:length(f_e),f_e-f);
E=mean(f_e)
MSE=mean((f_e-f).^2)