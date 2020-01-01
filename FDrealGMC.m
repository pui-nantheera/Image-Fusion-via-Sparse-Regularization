%Code of Joint image fusion and deconvolution on the real OCT and fundus
%image data sets.
clear all
% close all
% clc
getd = @(p)path(p,path);
getd('Metrics/');
addpath('utilities/');

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
    z1 = imresize(z1, [256 256]*2);
    z2 = imresize(z2, [256 256]*2);
end
z2 = z2(:,:,1);
z1 = z1(:,:,1);
A = uint8(z1*255);
B = uint8(z2*255);


curreseultDir = [resultDir, seq,'/'];
mkdir(curreseultDir);

% %read test image
% z1 = im2double(imread('./original image data/1_1_f.tif'));
% z1=z1(1:600,1:600);
%
% z2 = imread('./original image data/1_1_o.tif');
% z2= rgb2gray(z2);
% z2=im2double(z2);
% z2=z2(1:600,1:600);

%PSF estimation
psf_size=[5 5];
INITPSF = ones(psf_size);
[H1, PSF1] = estimatePSF(z1, psf_size);
for max_value=0.055
    PSF1=PSF1/sum(PSF1(:));%(PSF1/max(max(PSF1)))*max_value;
    [H2, PSF2] = estimatePSF(z2, psf_size);
    PSF2=PSF2/sum(PSF2(:));%(PSF2/max(max(PSF2)))*max_value;
%     [J, PSF1] = deconvblind(z1,INITPSF);
%     [J, PSF2] = deconvblind(z2,INITPSF);
    
    %Downsampling and upsampling
    K=1;
    [G1,Gt1]    = GGt(PSF1,K);
    [G2,Gt2]    = GGt(PSF2,K);
    
    %Wavelet and inverse wavelet transform
    levels = 4;
    [a,L] = wavedec2(z1,levels,'Haar');
    
    Wt  = @(x) wavedec2(x,levels,'Haar');
    Wtt = @(x) waverec2(x,L,'Haar');
    
    %Two sensor gain estimation methods
    [recovbete1, recovbete2] = SensorGain(z1, z2);
    % [recovbete1 recovbete2] = correctsensor(z1, z2);
    
    %Parameter setting
    for lam=0.0005%[0.0001 0.0005 0.001]
        rho=1;
        for gam=0.8%[0.6:0.1:1.3]
            tic
            miu=1.9 / ( rho * max( 1,  gam / (1-gam) ) );
            alp=lam*miu;
            [q,w]=size(z1);
            thk=zeros(1,q*w);
            v=zeros(1,q*w);
            itrs=20;
            for j=1:itrs
                
                w=thk-miu*(Wt(Gt1(recovbete2.*(recovbete2.*G1(Wtt(thk))-z1)))+Wt(Gt2(recovbete1.*(recovbete1.*G2(Wtt(thk))-z2)))+gam*(Wt(Gt1(recovbete2.*(recovbete2.*G1(Wtt(v-thk)))))+Wt(Gt2(recovbete1.*(recovbete1.*G2(Wtt(v-thk)))))));
                u=v-miu*gam*(Wt(Gt1(recovbete2.*(recovbete2.*G1(Wtt(v-thk)))))+Wt(Gt2(recovbete1.*(recovbete1.*G2(Wtt(v-thk))))));
                thk=(max(abs(w)-alp,0)).*sign(w);
                v=(max(abs(u)-alp,0)).*sign(u);
            end
            final=Wtt(thk);
            final = final/max(final(:));%*(0.5*(max(z1(:))+max(z2(:))));
            toc
            %Results showing
            peval = imqmet(cat(3, z1, z2), final);
            cvval = Cvejic_metric(cat(3, z1, z2),final);
            pvtrovic = petmetric(cat(3, z1, z2), final);
            imwrite(final,[curreseultDir, seq, 'real_','am',num2str(0),'.png']);
            disp(num2str([peval cvval pvtrovic]))
            figure
            imshow([z1 z2 final]); title([num2str(lam),', score=', num2str(pvtrovic),', ', num2str(peval),', ', num2str(cvval)])

        end
        
    end
end
