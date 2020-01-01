function [bete1, bete2] = SensorGain(VI, IR, p)

Sen_ir_inial = IR;
Sen_vi_inial = VI;

Sen_ir=Sen_ir_inial ;
Sen_vi=Sen_vi_inial ;

if nargin < 3
    p = 1; % block size
end
i = 90;
mbSize = p;

bete1 = zeros(size(Sen_ir));
bete2 = bete1;
row = size(Sen_ir,1);
col = size(Sen_ir,2);

for i = 1 : mbSize+1 : row-mbSize
    for j = 1 : mbSize+1 : col-mbSize
        
        
        R_ir = Sen_ir(i:i+p,j:j+p);
        R_vi = Sen_vi(i:i+p,j:j+p);
        
        f_ir =R_ir(:);%50* ones(size(R_ir(:)));
        
        f_vi =R_vi(:);%50* ones(size(R_vi(:)));
        
        vK = [f_ir,f_vi];
        variance = zeros(2); % 2 is the sensor numbers
        for k = 1:size(vK,1)
            
            temp_variance =  double(vK(k,:)') ;
            variance = variance+temp_variance*temp_variance';
            
        end
        variance;
        [s v]= eig(variance);
        
        % local_bete1=abs(v(1,1))
        % local_bete2=abs(v(2,2))
        local_bete1=abs(s(1,2));
        local_bete2=abs(s(2,2));
        if local_bete1 == local_bete2
            local_bete1 = 1;
            local_bete2 = 1;
        end
        bete1(i:i+mbSize,j:j+mbSize) = local_bete1.*ones(mbSize+1,mbSize+1);
        bete2(i:i+mbSize,j:j+mbSize) = local_bete2.*ones(mbSize+1,mbSize+1);
    end
end
%% normalizing

norm = bete1.^2 + bete2.^2;
%norm(find(norm==0))=1;

bete1 = bete1./(norm.^0.5); %% para for ir
bete2 = bete2./(norm.^0.5); %% para for vi

% bete1 = bete1./(bete1+bete2);
% bete2 = bete2./(bete1+bete2);

bete1(isnan(bete1))=1;
bete2(isnan(bete2))=1;




temp_fused = (double(Sen_ir)+double(Sen_vi))./2;

%figure,imshow(uint8(temp_fused))

% figure,imshow(uint8(temp_fused.*bete1))  % ir
% figure,imshow(uint8(temp_fused.*bete2))  % vi