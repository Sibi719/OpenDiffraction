clc;
clear all;
close all;
format long 
wavelength=500*10^-9;
m=2;
z1=0;
z2=0.5;
zi=z2/(m+1);
n=2^11;
R=1*10^-3;
delz2=z2-zi;
dimension=6*10^-3;
d1=dimension/n;
di=(wavelength*zi)/(n*d1);
d2=(wavelength*delz2)/(n*di);
F=(R^2)/(wavelength*z2);
disp(['Fresnel Number = ',num2str(F)]);
k=(2*pi/wavelength);
x1= (-n/2:((n/2)-1))*d1;
y1= (-n/2:((n/2)-1))*d1;
[X1,Y1]= meshgrid(x1,y1);
z_min=R*d1/wavelength
M=zeros(n);
A=(((X1.^2)+(Y1.^2))<=((R)^2));
M(A)=1;

figure
imagesc(x1*10^3,y1*10^3,M);
xlabel("x(mm)");
ylabel("y(mm)");
title("Source plane");
colormap(gray)
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
saveas(gcf,"aperture.png");

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
title("Observation plane")
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
saveas(gcf,"TS_z_"+z2+"_m_"+m+"_F_"+F+".png");
