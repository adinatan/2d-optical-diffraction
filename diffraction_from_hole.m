clear all
close all
i1=sqrt(-1);
range=24; %all distance units in millimeters
N=2^9;
 
lambda=2e-3;

 
 scale=linspace(-range/2,range/2,N);
ftscale=(N/range^2)*scale;
[x,y]=meshgrid(scale,scale);
z=linspace(0,1000,50); % propagation grid

[fx,fy]=meshgrid(ftscale,ftscale);
 %beam parametres
sigma=5;
 
A=exp(-(x.^2+y.^2)/(2*sigma^2));
hole_radius=2;
A_filtered=A.*(sqrt(x.^2+y.^2)>hole_radius); %beam after iris 
%A_filtered=A.*(sqrt(x.^2+y.^2)>0.9*iris_radius & %sqrt(x.^2+y.^2)<1.1*iris_radius); %ring beam for fun

%%
Aw=fftshift(fft2((A_filtered)));
z=linspace(0,5000,100);
dg=1/lambda^2-fx.^2-fy.^2;
k0=2*pi/lambda;
clims=[0 2]; % setting colormap min and max values
   
for j=1:length(z) ;
    if dg>0
      %  ff = exp(-i1*2*pi*z(j).*sqrt(dg)); %calculates the free space transfer for diffraction, There is no real advantage for the Fresnel approximation in numerical work, so we use the exact scalar transfer function.
        ff = exp(-i1*2*pi*(z(j).* sqrt(dg)-z(j)/lambda));
        % for large values of z(j), the large phase of this function tends to disrupt Fourier transforms.
        % substituting ff = exp(-i1*2*pi*(z(j).* sqrt(dg)-z(j)/lamda)); seems to correct this problem.
       % ff = exp(i1*(fx.^2+fy.^2)/(2*k0)*z(j));% fresnel approx (up to some factor , (2*pi?)).
    else
        ff=0;
    end
    
    ft=ff.*Aw;
    Eout=ifft2(fftshift(ft));
     h=imagesc(scale,scale,abs(Eout).^2,clims); hold on
    axis('square'); xlabel('mm') ; ylabel('mm')
    getframe;    
   
end
 
