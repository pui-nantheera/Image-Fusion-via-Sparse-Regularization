function [G,Gt] = GGt(h,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operators for super-resolution
% Stanley Chan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G  = @(x) fdown(x,h,K);
Gt = @(x) upf(x,h,K);
end

function y = fdown(x,h,K)
y= imfilter(x,h,'circular');

end

function y = upf(x,h,K)

y = imfilter(x,h','circular');

end
