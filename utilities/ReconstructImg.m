function [constructedImg, weightMtx] = ReconstructImg(imgPatchs, distance, a, b, n)
weightMtx = zeros(a,b);
constructedImg = zeros(a,b);
count = 0;
for i=1:distance:a-n+1
    for j = 1:distance:b-n+1
      count = count+1;
      imgPatch = reshape(imgPatchs(:,count),[n,n]); 
      constructedImg(i:i+n-1,j:j+n-1) = constructedImg(i:i+n-1,j:j+n-1) + imgPatch; 
      weightMtx(i:i+n-1,j:j+n-1) = weightMtx(i:i+n-1,j:j+n-1)+1; 
    end
end;