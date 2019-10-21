
clear all
close all
i1=sqrt(-1);
range=25; %all distance units in millimeters
npoints=2^2;

step=range/npoints;
lambda=800e-6;

scale=linspace(-range/2,range/2,npoints);
%x=scale; y=scale;
[x,y]=meshgrid(scale,scale);
%ftscale=(npoints/range^2)*scale;
%[fx,fy]=meshgrid(ftscale,ftscale);
sigma=10;
A=exp(-pi*(x.^2+y.^2)/(2*sigma^2));
iris_radius=3;
Ein=A.*(sqrt(x.^2+y.^2)<iris_radius);
for z=1
Eout=FresnelPropagator2(Ein,x,y,z,lambda,-1);
imagesc(abs(Eout).^2);
end
