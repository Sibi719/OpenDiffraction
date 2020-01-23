clc;
clear all;
close all;
Wavelength= 500*10^-9;
N=2^11;
k=(2*pi)/Wavelength;
R=1*10^-3;
Dimesnion=6*10^-3;
del1=Dimesnion/N;
z_min=R*del1/Wavelength
z=0.5;
del2=(Wavelength*z)/(N*del1)
M=zeros(N);
x1=((-N/2):((N/2)-1))*del1;
y1=((-N/2):((N/2)-1))*del1;
[X1,Y1]=meshgrid(x1,y1);
A=X1.^2+ Y1.^2<=R.^2;
M(A)=1;
F=R^2/(z*Wavelength);

figure
imagesc(x1*10^3,y1*10^3,M);
xlabel("x(mm)");
ylabel("y(mm)");
axis image
title(['Circular aperture of radius ',num2str(R*10^3),' mm']);
colormap(gray)
colorbar

x2=((-N/2):((N/2)-1))*del2;
y2=((-N/2):((N/2)-1))*del2;
df_1= 1/(N*del1);
exp_term= exp(((1j*k)/(2*z)).* (X1.^2+Y1.^2));
fft_in= M.*exp_term;
fft_ou= fftshift(fft2(fft_in));
[X2,Y2]=meshgrid(x2,y2);
expo=exp(((1j*k)/(2*z)).* (X2.^2+Y2.^2)) ;
irre=exp(1j.*k.*z)/(1j*Wavelength*z);
Output=irre*expo.*fft_ou;
scal=(sum(sum((abs(fft_in).^2)*(del1)*(del1))))/(sum(sum((abs(fft2(fft_in)).^2)*(df_1)*(df_1))));
I=abs(Output).^2;

figure
imagesc(x2*10^3,y2*10^3,I.*scal);
xlabel("x(mm)");
ylabel("y(mm)");
title("Observation plane");
colormap(hot)
colorbar
axis on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"OS_z_"+z+"_F_"+F+".png");

