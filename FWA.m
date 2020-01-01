%Code of image fusion with GMC regularization.
getd = @(p)path(p,path);
getd('Metrics/');
addpath('utilities/');

% %original image data
% z1 = im2double(imread('C:\Dropbox\Others\Rencheng Report\32951\Rencheng Zheng report and MATLAB source code\source code\original image data\ct.jpg'));
% z2 = im2double(imread('C:\Dropbox\Others\Rencheng Report\32951\Rencheng Zheng report and MATLAB source code\source code\original image data\mri.jpg'));
% z2 = z2(:,:,1);
% z1 = z1(:,:,1);

resultDir = '/Users/eexna/Work/Ultrasound imaging/Fusion/Rencheng Report/results/';
seq = 'slika';

%original image data
z2colour = [];
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
    z1 = imresize(z1, [1 1]*512);
    z2 = imresize(z2, [1 1]*512);
end
z2 = z2(:,:,1);
z1 = z1(:,:,1);
A = uint8(z1*255);
B = uint8(z2*255);

tic
[c1,s] =  wavedec2(z1,3,'Haar');
[c2,s] =  wavedec2(z2,3,'Haar');
final =  waverec2(0.5*(c1+c2),s,'Haar');
toc
peval = imqmet(cat(3, A, B), uint8(255*final));
cvval = Cvejic_metric(cat(3, A, B),uint8(255*final));
pvtrovic = petmetric(cat(3, A, B), uint8(255*final));
disp(num2str([peval cvval pvtrovic]))

figure
imshow(final); title(['score=', num2str(pvtrovic),', ', num2str(peval),', ', num2str(cvval)])
