%Code of image fusion with L1 regularization.
getd = @(p)path(p,path);
getd('Metrics/');
addpath('./utilities');

resultDir = '/Users/eexna/Work/Ultrasound imaging/Fusion/Rencheng Report/results/';
seq = 'slika';

%original image data
if strcmp(seq, 'clocks')
    z1 = im2double(imread('/Users/eexna/Work/Ultrasound imaging/Fusion/Rencheng Report/data/clocks/clockA.jpg'));
    z2 = im2double(imread('/Users/eexna/Work/Ultrasound imaging/Fusion/Rencheng Report/data/clocks/clockB.jpg'));
elseif strcmp(seq,'Bild2')
    commondir = '/Users/eexna/Work/Ultrasound imaging/Fusion/data/bild/';
    z1 = im2double(imread([commondir,'Bild2_ir_wh.bmp']));
    z2 = im2double(imread([commondir,'Bild2_tv.bmp']));
    z2 = rgb2ycbcr(z2);
    z2colour = z2;
elseif strcmp(seq, 'quad')
    commondir = '/Users/eexna/Work/Ultrasound imaging/Fusion/data/junction/';
    z1 = im2double(imread([commondir,'quad_TV_bw_user7_im3.jpg']));
    z2 = im2double(imread([commondir,'quad_IR_user4_im5.jpg']));
elseif  strcmp(seq, 'Octec')
    commondir = '/Users/eexna/Work/Ultrasound imaging/Fusion/data/Octec/';
    z1 = im2double(imread([commondir,'Octec_IR_22.jpg']));
    z2 = im2double(imread([commondir,'Octec_TV_22.jpg']));
    z2 = rgb2ycbcr(z2);
    z2colour = z2;
elseif strcmp(seq, 'slika')
    z1 = im2double(imread('/Users/eexna/Downloads/TestingImageDataset/testna_slika2a.bmp'));
    z2 = im2double(imread('/Users/eexna/Downloads/TestingImageDataset/testna_slika2b.bmp'));
    z1 = imresize(z1, [256 256]*2);
    z2 = imresize(z2, [256 256]*2);
end
z2 = z2(:,:,1);
z1 = z1(:,:,1);
A = uint8(z1*255);
B = uint8(z2*255);


% z1 = im2double(imread('./original image data/ct.jpg'));
% z1=[z1,zeros(159,1)];
% z1=[z1;zeros(1,160)];
% z2 = im2double(imread('./original image data/mri.jpg'));
% z2=[z2,zeros(159,1)];
% z2=[z2;zeros(1,160)];

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

% figure
% imshow(z1);
% figure
% imshow(z2);

% Two sensor gain estimation methods
% [recovbete1 recovbete2] = correctsensor(z1, z2);
[recovbete1, recovbete2] = SensorGain(z1, z2);


% Wavelet and inverse wavelet transform
levels = 4;
[a,L] = wavedec2(z1,levels,'Haar');

Wt  = @(x) wavedec2(x,levels,'Haar');
Wtt = @(x) waverec2(x,L,'Haar');
tic
%parameter setting
lam=0.001;
miu=1.9;
alp=lam*miu;
[q,w]=size(z1);
thk=zeros(1,size(a,2));
%main L1 regularization algorithm
itrs=25;
for j=1:itrs
    
    ck=thk-miu*(Wt(recovbete2.*(recovbete2.*Wtt(thk)-z1))+Wt(recovbete1.*(recovbete1.*Wtt(thk)-z2)));
    thk=(max(abs(ck)-alp,0)).*sign(ck);
    
end
final=Wtt(thk);
final = final/max(final(:))*(0.5*(max(z1(:))+max(z2(:))));
final_old = final;
toc


