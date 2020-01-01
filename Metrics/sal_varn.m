function s= sal_varn(img,block_size)
%
% sal_var.m Compute variance as Wang and Bovick
% 
%Input : 
%           img        : images 
%           block_siz  : (optional, default 8) size of the window 
%
%Output:    
%           s:  a variance saliency map of img. 
%
%========================================================================

if (nargin == 0 | nargin > 2)
   sl = -1*ones(size(img1));
   return;
end

if (nargin == 1)
   block_size = 8;
end

N = block_size.^2;
sum2_filter = ones(block_size);

img_sq   = img.*img;

img_sum    = (convn(img, sum2_filter, 'valid'))/N;
img_sq_sum = (convn(img_sq, sum2_filter, 'valid'))/N;

s= img_sq_sum - (img_sum.*img_sum);










