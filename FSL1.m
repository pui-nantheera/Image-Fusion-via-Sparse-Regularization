%Code of joint image fusion and super-resolution with L1 regularization.
clear all
close all
clc
getd = @(p)path(p,path); % scilab users must *not* execute this
getd('Metrics/');
addpath(genpath('./utilities/'));


%original image data
z1 = im2double(imread('./original image data/ct.jpg'));
z1=[z1,zeros(159,1)];
z1=[z1;zeros(1,160)];
z2 = im2double(imread('./original image data/mri.jpg'));
z2=[z2,zeros(159,1)];
z2=[z2;zeros(1,160)];

% z1 = imread('./original image data/011_CT_1.bmp');
% z2 = imread('./original image data/011_T2.bmp');
% z1= rgb2gray(z1);
% z2= rgb2gray(z2);
% z1=im2double(z1);
% z2=im2double(z2);

% z1 = im2double(imread('./original image data/1_1_f.tif'));
% z1=z1(1:600,1:600);
% z2 = imread('./original image data/1_1_o.tif');
% z2= rgb2gray(z2);
% z2=im2double(z2);
% z2=z2(1:600,1:600);

%Point spread function
h = fspecial('gaussian',[9 9],1);
%Scale parameter
K = 2;
%calculate the observed image
y1 = imfilter(z1,h,'circular');
y1 = downsample2(y1,K);
y2 = imfilter(z2,h,'circular');
y2 = downsample2(y2,K);

%Super-resolution operator
[rows_in,cols_in] = size(y1);
rows      = rows_in*K;
cols      = cols_in*K;
N         = rows*cols;
[G,Gt]    = defGGt(h,K);
GGt       = constructGGt(h,K,rows,cols);

%Wavelet and inverse wavelet transform
L=[40,40;40,40;80,80;160,160];
% L=[64,64;64,64;128,128;256,256];
% L=[150,150;150,150;300,300;600,600];
Wt  = @(x) wavedec2(x,2,'Haar');
Wtt = @(x) waverec2(x,L,'Haar');

%Two sensor gain estimation methods
[recovbete1 recovbete2] = SensorGain(y1, y2);
% [recovbete1 recovbete2] = correctsensor(y1, y2);


lam=0.0005;
miu=1.9;
alp=lam*miu;
[q,w]=size(z1);
thk=zeros(1,q*w);

%main L1regularization algorithm
itrs=50;
for j=1:itrs
    
    ck=thk-miu*(Wt(Gt(recovbete2.*(recovbete2.*G(Wtt(thk))-y1)))+Wt(Gt(recovbete1.*(recovbete1.*G(Wtt(thk))-y2))));
    %     ck=thk-miu*((thk-noisyPatchs1(:, i))+(thk-noisyPatchs2(:, i)));
    thk=(max(abs(ck)-alp,0)).*sign(ck);
end
final=Wtt(thk);

%Results showing
figure
imshow(final);
