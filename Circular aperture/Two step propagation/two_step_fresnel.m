clc;
clear all;
format long 
n=2^12;
R=1*10^-3;
wavelength=500*10^-9;
dimension=10*10^-3;
d1=dimension/n;
zmin=n*d1*d1/wavelength
m=2;
z1=0;
z2=0.5;
zi=z2/(m+1);
delz2=z2-zi;
di=(wavelength*zi)/(n*d1);
d2=(wavelength*delz2)/(n*di);
k=(2*pi/wavelength);
x1= (-n/2:((n/2)-1))*d1;
y1= (-n/2:((n/2)-1))*d1;
[X1,Y1]= meshgrid(x1,y1);
M=zeros(n);
A=(((X1.^2)+(Y1.^2))<=((R)^2));
M(A)=1;
figure
imagesc(x1*10^3,y1*10^3,M);
axis image
title(['Single circular aperture of radius ',num2str(R),'m']);
colorbar

xi= ((-n/2):((n/2)-1))*di;
yi= ((-n/2):((n/2)-1))*di;
[Xi,Yi]=meshgrid(xi,yi);

exp_term=exp(((1j*k)/(2*zi)).*(X1.^2+Y1.^2));
fft_in= M.*exp_term;
fft_ou= fftshift(fft2(fft_in));
expo=exp(((1j*k)/(2*zi)).*(Xi.^2+Yi.^2)) ;
irre=exp(1j.*k.*zi)/(1j*wavelength*zi);
df_1= 1/(n*d1);
scal=sqrt((sum(sum((abs(fft_in).^2)*(d1)*(d1))))/(sum(sum((abs(fft2(fft_in)).^2)*(df_1)*(df_1)))));
Outputi=irre*expo.*fft_ou.*scal;

x2= (-n/2:((n/2)-1))*d2;
y2= (-n/2:((n/2)-1))*d2;
[X2,Y2]= meshgrid(x2,y2);
exp_term=exp(((1j*k)/(2*delz2)).*(Xi.^2+Yi.^2));
fft_in= Outputi.*exp_term;
fft_ou= (fft2(fft_in));
expo=exp(((1j*k)/(2*delz2)).* (X2.^2+Y2.^2)) ;
irre=exp(1j.*k.*delz2)/(1j*wavelength*delz2);
Output=irre*expo.*fft_ou;
df_i= 1/(n*di);
scal=((sum(sum((abs(fft_in).^2)*(di)*(di))))/(sum(sum((abs(fft2(fft_in)).^2)*(df_i)*(df_i)))));
I=abs(Output).^2.*scal;

figure
imagesc(x2*10^3,y2*10^3,I);
xlabel("x(mm)");
ylabel("y(mm)");
colormap(hot)
colorbar
