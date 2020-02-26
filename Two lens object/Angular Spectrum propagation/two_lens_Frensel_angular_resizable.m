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
m=2;
del2=m*del1;
zmax=N*del1*del2/Wavelength
z=zmax*2;
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

x2=((-N/2):((N/2)-1))*del2;
y2=((-N/2):((N/2)-1))*del2;
[X2,Y2]=meshgrid(x2,y2);

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

df_i=1/(N*del1);
fx1=((-N/2):((N/2)-1))*df_i;
fy1=((-N/2):((N/2)-1))*df_i;
[Fx1,FY1]=meshgrid(fx1,fy1);
A=Fx1.^2+ FY1.^2;
%phs= exp(1j*k*z).*exp(1j*k*(X1.^2+Y1.^2)/(2*z))./(1j*z*Wavelength);
phs= exp(-1j*pi*z*Wavelength*(A)/m);
F1=fftshift(fft2((M.*exp((1j.*k.*(1-m).*(X1.^2+Y1.^2))./(2.*z)))/m));
F2=fftshift(fft2(phs));
ph2=exp((-1j.*k.*((1-m)/m).*(X2.^2+Y2.^2))./(2.*z));

UU= ph2.*(ifft2(ifftshift(F1.*phs)));
I=(abs(UU).^2);

figure
imagesc(x2*10^3,y2*10^3,I);
xlabel("x(mm)");
ylabel("y(mm)");
colorbar
colormap(jet)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Papers\Open_Diffraction\ang1.png");