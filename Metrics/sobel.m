function [G,A]=sobel(imagein,regi)

if nargin < 2
    regi='same';
elseif ~ (strcmp(regi,'valid') | strcmp(regi,'same'))
    error('Wrong input argument')
end

[r c n] = size(imagein);

h1=[[1 2 1];[0 0 0];[-1 -2 -1]]/8;
h3=[[-1 0 1];[-2 0 2];[-1 0 1]]/8;

for i = 1:n
    x(:,:,i)=filter2(h1,imagein(:,:,i),regi);
    y(:,:,i)=filter2(h3,imagein(:,:,i),regi);
end
G=(x.^2+y.^2).^(.5);
A=atan(y./(x+eps));
A(isnan(A))=0;