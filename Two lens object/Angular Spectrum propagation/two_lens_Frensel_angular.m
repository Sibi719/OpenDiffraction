clc;
clear all;
close all;
Wavelength=500*10^-9;
N=2^12;
k=(2*pi)/Wavelength;
R1=2.3*10^-2;
R2=1.4*10^-2;
Dimesnion=8*10^-2;
del1=Dimesnion/N;
z=2;
zmax=N*del1*del1/Wavelength
alpha1=0.0001;
beta1= 0.000002;
alpha2=alpha1;
beta2=beta1;
cx1=0;
cy1=0;
cx2=-1.5*10^-2;
cy2=0;
M=zeros(N);
x1=((-N/2):((N/2)-1))*del1;
y1=((-N/2):((N/2)-1))*del1;
[X1,Y1]=meshgrid(x1,y1);

L1= real(sqrt(R1^2-((X1-cx1).^2 + (Y1-cy1).^2)));
L2= real(sqrt(R2^2-((X1-cx2).^2 + (Y1-cy2).^2)));

phase1 = k*alpha1*L1;
phase2 = k*alpha2*L2;
phase= phase1+phase2;

T1 = k*beta1*L1;
T2 = k*beta2*L2;
T= (exp(-(T1+T2))).^2;

figure
imagesc(x1*10^3,y1*10^3,T);
axis image
xlabel("x(mm)");
ylabel("y(mm)");
colorbar
colormap(jet)

figure
imagesc(x1*10^3,y1*10^3,phase);
axis image
xlabel("x(mm)");
ylabel("y(mm)");
colorbar
colormap(jet)

M=sqrt(T).*exp(1j.*phase);

[IP]=FresProp(del1,z,Wavelength,N,M);

figure
imagesc(x1*10^3,y1*10^3,IP);
axis image
xlabel("x(mm)");
ylabel("y(mm)");
colorbar
colormap(jet)

function [I]=FresProp(dpix,z,lambda,Hsize,Hcrop)
 
%Spatial frequencies
Xsize = Hsize*dpix; %Hsize is the number of pixel, dpix is the length of one pixel, Xsize is the total lenght of the image. 
du = 1/(Xsize);% the resolution of fourier frequency coordinates
%Nyquist cut-off for Sampling Hologram
umax = 1/(2*dpix); %define the k space 
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
clear  u V  du;
 
%Evanescent cut-off 
uev = 1/lambda; %???????
 
%Nyquist cut-off for Fresnel Propagation Kernel
unp = uev*(Xsize/(2*abs(z)));
clear Xsize;
 
%Circular window
A = U.^2+(U').^2;
clear U;
if uev>=unp
    ucut = unp;
end
if unp>uev
    ucut = uev;
end
W= sqrt(A);
W = (W<=ucut); 
% disp(['Cutoff =',num2str(ucut),' Evansecent Cutoff =',num2str(uev),...
%' Nyquist Cutoff =', num2str(unp),'u max =',num2str(umax)])
clear ucut uev unp
 
%Fresnel kernel: paraxial approximation
H = exp((-1i*pi*lambda* z).*(A));
clear A;
 
%Truncate kernel
H = W.*H;
clear W;
 
%Hologram Spectrum
Htemp = fft2(Hcrop);
HH = fftshift(Htemp);
clear Htemp;
 
%Propagate field
RR = HH.*H;
clear H HH;
RR =ifftshift(RR);
R = ifft2(RR);
I=abs(R).^2;
clear RR;
end