function [SF] = spfreqmat (mat,block,shape)
%
%spfreq.m,v 1.11 2002/03/04 17:37:40
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2003
%===========================================================================
%
%       	[SF] = spfreqmat (mat,block,shape)              
%
%	Calculates the spacial frequency of the input matrix
%
%===========================================================================
%
%   
if nargin < 3
    shape='valid';
elseif ~ (strcmp(shape,'valid') | strcmp(shape,'same'))
    error('Wrong input argument')
end
    
rff=[1,-1];
RF=(convn(mat,rff,shape)).^2;
cff=[1;-1];
CF=(convn(mat,cff,shape)).^2;

NRF=ones(block,block-1);
RFS=convn(RF,NRF,shape);
NFS=ones(block-1,block);
CFS=convn(CF,NFS,shape);

SF=sqrt(RFS+CFS)/block;