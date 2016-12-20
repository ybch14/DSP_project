close all;clear ;clc;
load music.mat;
f_music=f_e;
load esprit.mat;
f_esprit=f_e;
figure;
subplot(2,1,1)
plot(1:length(f_music),f_music-324.2333);
title('MUSIC algorithm error');
xlabel('time');
ylabel('Hz^2');
subplot(2,1,2)
plot(1:length(f_esprit),f_esprit-324.2333);
title('ESPRIT algorithm error');
xlabel('time');
ylabel('Hz^2');

