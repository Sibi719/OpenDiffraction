clc;
clear all;
close all;
Wavelength= 500*10^-9;
z=0.03;
N=2^12;
k=(2*pi)/Wavelength;
R=1*10^-3;
Dimesnion=10*10^-3;
del1=Dimesnion/N;
del2=(Wavelength*z)/(N*del1);
M=zeros(N);
x1=((-N/2):((N/2)-1))*del1;
y1=((-N/2):((N/2)-1))*del1;
m=1;
del2=m*del1;
x2=((-N/2):((N/2)-1))*del2;
y2=((-N/2):((N/2)-1))*del2;
[X1,Y1]=meshgrid(x1,y1);
[X2,Y2]=meshgrid(x2,y2);
A=X1.^2+ Y1.^2<=R.^2;
M(A)=1;
F=R^2/(z*Wavelength)
%Aperture
figure
imagesc(x1*10^3,y1*10^3,M);
xlabel("x(mm)");
ylabel("y(mm)");
axis image
colorbar
colormap(gray)


df_i=1/(N*del1);
fx1=((-N/2):((N/2)-1))*df_i;
fy1=((-N/2):((N/2)-1))*df_i;
[Fx1,FY1]=meshgrid(fx1,fy1);
A=Fx1.^2+ FY1.^2;
phs= exp(-1j*pi*z*Wavelength*(A)/m);
F1=fftshift(fft2((M.*exp((1j.*k.*(1-m).*(X1.^2+Y1.^2))./(2.*z)))/m));
F2=fftshift(fft2(phs));
ph2=exp((-1j.*k.*((1-m)/m).*(X2.^2+Y2.^2))./(2.*z));

UU=exp(1j.*k.*z).* ph2.*(ifft2(ifftshift(F1.*phs)));
I=(abs(UU).^2);

figure
imagesc(x2*10^3,y2*10^3,I);
xlabel("x(mm)");
ylabel("y(mm)");
colorbar
colormap(hot)
