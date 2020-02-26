clc;
clear all;
close all;
format long 
n=2^12;
wavelength=500*10^-9;
dimension=10*10^-3;
d1=dimension/n;
zmin=n*d1*d1/wavelength
m=1;
z1=0;
z2=zmin;
zi=z2/(m+1);
delz2=z2-zi;
di=(wavelength*zi)/(n*d1);
d2=(wavelength*delz2)/(n*di);
R1=2.3*10^-3;
R2=1.4*10^-3;
k=(2*pi/wavelength);
alpha1=0.001;
beta1=0.00002;
alpha2=alpha1;
beta2=beta1;
cx1=0;
cy1=0;
cx2=-1.5*10^-3;
cy2=0;
x1= (-n/2:((n/2)-1))*d1;
y1= (-n/2:((n/2)-1))*d1;
[X1,Y1]= meshgrid(x1,y1);
M=zeros(n);
L1= real(sqrt(R1^2-((X1-cx1).^2 + (Y1-cy1).^2)));
L2= real(sqrt(R2^2-((X1-cx2).^2 + (Y1-cy2).^2)));

phase1 = k*alpha1*L1;
phase2 = k*alpha2*L2;

phase= phase1+phase2;

T1 = k*beta1*L1;
T2 = k*beta2*L2;

T= (exp(-(T1+T2))).^2;
% 
% figure
% imagesc(x1*10^3,y1*10^3,T);
% axis image
% xlabel("x(mm)");
% ylabel("y(mm)");
% colorbar
% colormap(jet)
% 
% figure
% imagesc(x1*10^3,y1*10^3,phase);
% axis image
% xlabel("x(mm)");
% ylabel("y(mm)");colorbar
% colormap(jet)

I_o=ones(n);

M=sqrt(I_o.*T).*exp(1j.*phase);

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
colormap(jet)
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"D:\Papers\Open_Diffraction\Figures\6.png");