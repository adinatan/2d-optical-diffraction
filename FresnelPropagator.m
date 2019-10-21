function Eout=FresnelPropagator(Ein,x,z,lambda,kcut)
%15/9/08
% Calculate the field E(x,z) for a given E(x,z=0)=Ein, according to Fersnel propagation. 
% lambda is the wavelength, and kcut is the cutoff for the spatial
% frequency 
% 
if nargin<5
    kcut=-1;
end

k0=2*pi/lambda;
dx=x(2)-x(1);
kx=linspace(-pi/dx,pi/dx,length(x));

Ein_k=fftshift(fft(Ein));
if kcut~=-1 
    Eout_k=exp(i*kx.^2/(2*k0)*z).*Ein_k.*exp(-(kx/kcut).^2);
else
    Eout_k=exp(i*kx.^2/(2*k0)*z).*Ein_k;
end
Eout=ifft(ifftshift(Eout_k));
