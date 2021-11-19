% ------- main procedure: call for LUTs to perform diffration 
% ------- calcualtion based on polygon-based method.
% ------- By Fan Wang, Chiba University, Japan, 2021/11/19.
% ------- wangfan@chiba-u.jp, phy_wangfan@163.com
tic
clear 
% clearvars -except  parameter LU2 LV2 LU0 LV0 LV1 LU1 
close all
%% with PCA
load('E:\hologramLUTs\coeff_LU0.mat');
load('E:\hologramLUTs\coeff_LU1.mat');
load('E:\hologramLUTs\coeff_LV0.mat');
load('E:\hologramLUTs\coeff_LV1.mat');
load('E:\hologramLUTs\score_LU0.mat');
load('E:\hologramLUTs\score_LU1.mat');
load('E:\hologramLUTs\score_LV0.mat');
load('E:\hologramLUTs\score_LV1.mat');
load('E:\hologramLUTs\meanLU0.mat');
load('E:\hologramLUTs\meanLU1.mat');
load('E:\hologramLUTs\meanLV0.mat');
load('E:\hologramLUTs\meanLV1.mat');
coeff_LU0=gpuArray(coeff_LU0);
score_LU0=gpuArray(score_LU0);
coeff_LV0=gpuArray(coeff_LV0);
score_LV0=gpuArray(score_LV0);
coeff_LV1=gpuArray(coeff_LV1);
score_LV1=gpuArray(score_LV1);
coeff_LU1=gpuArray(coeff_LU1);
score_LU1=gpuArray(score_LU1);
meanLU0=gpuArray(meanLU0);
meanLV0=gpuArray(meanLV0);
meanLU1=gpuArray(meanLU1);
meanLV1=gpuArray(meanLV1);
  
LU0=score_LU0*coeff_LU0'+meanLU0;
LV0=(score_LV0*coeff_LV0'+meanLV0)';
LU1=score_LU1*coeff_LU1'+meanLU1;
LV1=(score_LV1*coeff_LV1'+meanLV1)';

load('E:\hologramLUTs\LU2.mat');
load('E:\hologramLUTs\LV2.mat');
LU2=gpuArray(LU2);
LV2=gpuArray(LV2);
toc % elapsed time of loading data
load('E:\hologramLUTs\parameter_plg.mat');
%% import data
% t=clock;
tic
load('index_face6320_occlusion.mat'); %
load('vertex6320_occlusion.mat') %

Lx=parameter.Lx; 
Ly=parameter.Ly;
depth=parameter.depth;

max1=max(vertex); min1=min(vertex);
xmax=max1(1)-min1(1);ymax=max1(2)-min1(2);zmax=max1(3)-min1(3);
scaling_x=(Lx-1.5)/xmax; scaling_y=(Ly-1.5)/ymax; scaling_z=depth/zmax; % 缩放比例
vertex(:,1)=vertex(:,1)*scaling_x; 
vertex(:,2)=vertex(:,2)*scaling_y; 
vertex(:,3)=vertex(:,3)*scaling_z;
vertex(:,1)=vertex(:,1)-(max(vertex(:,1))+min(vertex(:,1)))/2;
vertex(:,2)=vertex(:,2)-((max(vertex(:,2))+min(vertex(:,2)))/2);
vertex(:,3)=vertex(:,3)-min(vertex(:,3));%((max(vertex(:,3))+min(vertex(:,3)))/2);
Lxo=max(vertex(:,1))-min(vertex(:,1));
Lyo=max(vertex(:,2))-min(vertex(:,2));
%% fundamental parameters
eps=parameter.eps;
lamda=parameter.lamda;
k=2*pi/lamda;
smp_x=parameter.smp_x;
smp_y=parameter.smp_y;
SMP_x=parameter.SMP_x;
SMP_y=parameter.SMP_y;
d=parameter.d;
prcs_xy=parameter.prcs_xy;
prcs_z=parameter.prcs_z;
dp_x=parameter.dp_x;
dp_y=parameter.dp_y;
d=ceil(d/prcs_z)*prcs_z;
vertex(:,1:2)=ceil(vertex(:,1:2)/prcs_xy)*prcs_xy;
vertex(:,3)=ceil(vertex(:,3)/prcs_z)*prcs_z+d;
% vertex(:,3)=ceil(vertex(:,3)/lamda)*lamda;


prcs_uv=parameter.prcs_uv;
dp_u0=round(1/Lx,abs(log10(prcs_uv)));   % precision of frequency u0
dp_v0=round(1/Ly,abs(log10(prcs_uv))); 
stp_u0=int32(dp_u0/prcs_uv);
stp_v0=int32(dp_v0/prcs_uv);
u0=-dp_u0/2-(smp_x/2-1)*dp_u0:dp_u0:dp_u0/2+(smp_x/2-1)*dp_u0;
v0=(dp_v0/2+(smp_y/2-1)*dp_v0:-dp_v0:-dp_v0/2-(smp_y/2-1)*dp_v0)';
[u0,v0]=meshgrid(u0,v0);
u0=gpuArray(u0);
v0=gpuArray(v0);
w0=1/lamda-lamda/2.*(u0.^2+v0.^2);

len_fx1=parameter.len_fx1;
len_fy1=parameter.len_fy1;
%% LUT method
FH=zeros(smp_y,smp_x);
for i=1:size(index_face,1)  %[29](397,1144)
    %% solving affine matrix T
    A=vertex(index_face(i,1),:);
    B=vertex(index_face(i,2),:);
    C=vertex(index_face(i,3),:);
    P=[0 1 1;0 0 1;0 0 0;1 1 1];  % corrdinates of primitive triangle
    V=[A' B' C';1 1 1]; % corrdinates of arbitrary triangle
    X=V(1,:);Y=V(2,:);Z=V(3,:);
    T=[X(2)-X(1) X(3)-X(2) 0 X(1);
       Y(2)-Y(1) Y(3)-Y(2) 0 Y(1);
       Z(2)-Z(1) Z(3)-Z(2) 0 Z(1)];
    du=T(3,1)/lamda;
    dv=T(3,2)/lamda;
   
    fu1_N=int16((T(1,1)+T(1,2)-floor(-Lx/10/prcs_xy)*prcs_xy)/prcs_xy+1); % exp(-j2π*(u'+v'))
    fv1_N=int16((T(2,1)+T(2,2)-floor(-Ly/10/prcs_xy)*prcs_xy)/prcs_xy+1);
    fz1_N=int16((T(3,1)+T(3,2)-floor(-depth/10/prcs_z)*prcs_z)/prcs_z+1);
    c1=CT(T(1,1)+T(1,2),T(2,1)+T(2,2),T(3,1)+T(3,2),lamda)*exp(1j*2*pi*(du+dv));
    
    fu2_N=int16((T(1,1)-floor(-Lx/10/prcs_xy)*prcs_xy)/prcs_xy+1); % exp(-j2π*(u'))
    fv2_N=int16((T(2,1)-floor(-Ly/10/prcs_xy)*prcs_xy)/prcs_xy+1);
    fz2_N=int16((T(3,1)-floor(-depth/10/prcs_z)*prcs_z)/prcs_z+1);
    c2=CT(T(1,1),T(2,1),T(3,1),lamda)*exp(1j*2*pi*(du));
          
    N_a34=int32((T(3,4)-d)/prcs_z+1);
    N_ua=int32(T(1,4)/T(3,4)/lamda/prcs_uv);
    N_vb=int32(T(2,4)/T(3,4)/lamda/prcs_uv);
    c0=exp(-1j*2*pi/lamda/2*(T(1,4)^2+T(2,4)^2)/T(3,4)); % CT(T(1,4),T(2,4),T(3,4),lamda);%
    J=abs(T(1,1)*T(2,2)-T(1,2)*T(2,1)); 
%% division term
    u=T(1,1).*u0+T(2,1).*v0+T(3,1).*w0-du;
    v=T(1,2).*u0+T(2,2).*v0+T(3,2).*w0-dv;
    D1=(-4*pi*pi).*(v.*(u+v));
    D2=(4*pi*pi).*(u.*v);
%% exponential term
    E1=LV1(:,(fz1_N-1)*len_fy1+fv1_N)*c1*LU1((fz1_N-1)*len_fx1+fu1_N,:)-1; % exp(-j2π*(u'+v'))-1
    E2=LV1(:,(fz2_N-1)*len_fy1+fv2_N)*c2*LU1((fz2_N-1)*len_fx1+fu2_N,:)-1; % exp(-j2π*(u'))-1

    if V(3,1)==V(3,3) || V(3,1)==V(3,2)
        if V(3,1)==V(3,3) && V(3,1)==V(3,2) 
            E1=LV2(:,fv1_N)*LU2(fu1_N,:)-1;
            E2=LV2(:,fv2_N)*LU2(fu2_N,:)-1;
        elseif  V(3,1)==V(3,3)  % when a31+a32=z3-z1=0  
            E1=LV2(:,fv1_N)*LU2(fu1_N,:)-1; % exp(-j2π*(u'+v'))                   
        elseif  V(3,1)==V(3,2) % when a31=z2-z1=0 
            E2=LV2(:,fv2_N)*LU2(fu2_N,:)-1; % exp(-j2π*(u'))      
        else
        end
    end    
    F=E1./D1+E2./D2;
    F = ifc(F,u,v); % if condition
    %% e^{-1j2π*[(a11+a12)u0+(a21+a22)v0+(a31+a32)w0]}
    lv0=LV0(size(LV0,1)/2-SMP_y/2+1+N_vb:stp_v0:size(LV0,1)/2+SMP_y/2+N_vb,N_a34);
    lu0=LU0(N_a34,size(LU0,2)/2-SMP_x/2+1-N_ua:stp_u0:size(LU0,2)/2+SMP_x/2-N_ua);
    E0=lv0*c0.*lu0;
    FH=FH+J*E0.*F;
end
toc % elapsed time of CGH calculation
holo=fftshift(ifft2(fftshift(FH))); % complex hologram
figure,imshow(abs(holo),[]),title('Holgoram');

di=d;
iF=FH.*(exp(1j*k*di*sqrt(1-(lamda*u0).^2-(lamda*v0).^2)));  %iH
iU=fftshift(ifft2(iF));
figure,imshow(abs(iU),[]),title('reconstruction by the complex holgoram');
figure,imshow(angle(iU),[]),title('phase reconstruction');

%% reference light
x=-dp_x/2-(smp_x/2-1)*dp_x:dp_x:dp_x/2+(smp_x/2-1)*dp_x;
y=dp_y/2+(smp_y/2-1)*dp_y:-dp_y:-dp_y/2-(smp_y/2-1)*dp_y;
[x,y]=meshgrid(x,y);
alpha_y=3; % off-axis along y dirction  % user-defined the angle of off-axis
alpha_x=0; % off-axis along y dirction
R=exp(1j*(k*sind(alpha_x)*x+k*sind(alpha_y)*y));
%% double phase-reconstruction D-fft
dpH=dphase(holo.*R); % double phase hologram interfered with reference light R
figure,imshow(dpH,[]),title('double phase coding hologram');
di=d+1;
iF=fft2(exp(1j*2*pi*dpH)).*exp(1j*k*di*(1-lamda^2/2*(u0.^2+v0.^2)));
Ii=ifft2(iF);
figure,imshow(abs(Ii),[]),title('reconstruction by the double phase coding hologram');

