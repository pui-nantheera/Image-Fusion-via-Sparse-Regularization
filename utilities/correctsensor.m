function [recovbete1,recovbete2]=correctsensor(noisyImg1,noisyImg2)
distance=1;
n=9;
[a,b] = size(noisyImg1);
noisyPatchs1 = zeros(n*n, ceil(1/distance*(a-n+1))*ceil(1/distance*(b-n+1)));
noisyPatchs2 = zeros(n*n, ceil(1/distance*(a-n+1))*ceil(1/distance*(b-n+1)));
count = 0;
for i = 1:distance:a-n+1
    for j = 1:distance:b-n+1
       count = count + 1;
       imgPatch1 = noisyImg1(i:i+n-1, j:j+n-1);
       noisyPatchs1(:, count) = imgPatch1(:);
    end
end
count = 0;
for i = 1:distance:a-n+1
    for j = 1:distance:b-n+1
       count = count + 1;
       imgPatch2 = noisyImg2(i:i+n-1, j:j+n-1);
       noisyPatchs2(:, count) = imgPatch2(:);
    end
end
for i=1:size(noisyPatchs1,2);
 vk=0;
    for j=1:81
        vk=vk+[noisyPatchs1(j,i),noisyPatchs2(j,i)]'*[noisyPatchs1(j,i),noisyPatchs2(j,i)];
    end
    vk=(1/80)*vk;
    [v,d] = eig(vk);
    value=max(d(:));
    [m n]=size(v);
    he=0;
    for p=1:m
    he=he+v(p,1);
    end
    vector=v(:,1);
    vector=abs(vector);
%     vector=abs(vector./norm(vector));
    beta(:,i)=vector;
    bete1(:,i)=beta(1,i).*ones(81,1);
    bete2(:,i)=beta(2,i).*ones(81,1);
end


n=9;
[recovbete1, weightMtx1] = ReconstructImg(bete1, distance, a, b, n);
recovbete1 = recovbete1./weightMtx1;
[recovbete2, weightMtx2] = ReconstructImg(bete2, distance, a, b, n);
recovbete2 = recovbete2./weightMtx2;

no=recovbete1.^2+recovbete2.^2;
recovbete1=recovbete1./(no.^0.5);
recovbete2=recovbete2./(no.^0.5);

end