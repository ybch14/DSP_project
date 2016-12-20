function [omega,ww,Pw]=my_music(y,n,m,delta_omega)
% y --- signal
% n --- number of sine signals
% m --- the size of convariance matrix
% delta_omega --- the resolution of frequency
N=length(y);
% compute the convariance matrix R
C=xcorr(y,'biased');
R=C(N:N+m-1)';
R=toeplitz(R);
% eigendecomposition, V is eigenvector 
% and D is diag(eigenvalues)
% the eigenvalues from small to big
[V,~]=eig(R);
% disp(D);
Nw=delta_omega;ww=linspace(0,2*pi,Nw+1);
e=exp(-1i*ww'*(0:m-1));
G=V(:,1:(m-n));
eg=e*G;
Pw=1./real(diag(eg*eg')');
[~,locs]=findpeaks(Pw,ww,'Threshold',1e-4);
omega=locs;
end