function Eout=FresnelPropagator2(Ein,x,y,z,lambda,kcut)
%15/9/08
% Calculate the field E(x,z) for a given E(x,z=0)=Ein, according to Fersnel propagation. 
% lambda is the wavelength, and kcut is the cutoff for the spatial
% frequency 
% 
if nargin<6
    kcut=-1;
end

k0=2*pi/lambda;
dx=x(2)-x(1);
kx=linspace(-pi/dx,pi/dx,length(x));

dy=y(2)-y(1);
ky=linspace(-pi/dy,pi/dy,length(y));

[KX,KY]=meshgrid(kx,ky);

Ein_k=fftshift(fft2(Ein));
if kcut~=-1 
    Eout_k=exp(i*(KX.^2+KY.^2)/(2*k0)*z).*Ein_k.*exp(-((KX.^2+KY.^2)/kcut^2));
else
    Eout_k=exp(i*(KX.^2+KY.^2)/(2*k0)*z).*Ein_k;
end
Eout=ifft2(ifftshift(Eout_k));
