clc;
clear all;
close all;
Wavelength= 500*10^-9;
N=2^12;
k=(2*pi)/Wavelength;
R=1*10^-3;
Dimesnion=10*10^-3;
del1=Dimesnion/N;
zmin= 2*R*del1/Wavelength
z=0.06;
del2=(Wavelength*z)/(N*del1);
M=zeros(N);
x1=((-N/2):((N/2)-1))*del1;
y1=((-N/2):((N/2)-1))*del1;
[X1,Y1]=meshgrid(x1,y1);
A=abs(X1)<=R & abs(Y1)<=R ;
M(A)=1;
figure
imagesc(x1*10^3,y1*10^3,M);
axis image

x2=((-N/2):((N/2)-1))*del2;
y2=((-N/2):((N/2)-1))*del2;
df_1= 1/(N*del1);
exp_term= exp(((1j*k)/(2*z)).* (X1.^2+Y1.^2));
fft_in= M.*exp_term;
fft_ou= fftshift(fft2(fft_in));

[X2,Y2]=meshgrid(x2,y2);
expo=exp(((1j*k)/(2*z)).* (X2.^2+Y2.^2)) ;
irre=exp(1j.*k.*z)/(1j*Wavelength*z);
scal=sqrt((sum(sum((abs(fft_in).^2)*(del1)*(del1))))/(sum(sum((abs(fft2(fft_in)).^2)*(df_1)*(df_1)))));
Output=irre*expo.*fft_ou.*scal;

I=abs(Output).^2;

figure
imagesc(x2*10^3,y2*10^3,I);
xlabel("x(mm)");
ylabel("y(mm)");
colormap(hot)
colorbar
