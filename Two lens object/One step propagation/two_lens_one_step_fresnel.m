clc;
clear all;
close all;
Wavelength=500*10^-9;
N=2^12;
k=(2*pi)/Wavelength;
R1=2.3*10^-3;
R2=1.4*10^-3;
Dimesnion=10*10^-3;
del1=Dimesnion/N;
zmin=N*del1*del1/Wavelength
z=0.1
del2=(Wavelength*z)/(N*del1);


alpha1=0.001;
beta1= 0.00002;
alpha2=alpha1;
beta2=beta1;
cx1=0;
cy1=0;
cx2=-1.5*10^-3;
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
xlabel("x(mm)");
ylabel("y(mm)");
colorbar
colormap(jet)



figure
imagesc(x1*10^3,y1*10^3,phase);
xlabel("x(mm)");
ylabel("y(mm)");
colorbar
colormap(jet)



M=sqrt(T).*exp(1j.*phase);

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
colormap(jet)
colorbar


