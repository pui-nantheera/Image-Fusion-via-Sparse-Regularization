function [sx,sy]=sobel_xy(a,varargin);

% Racunanje sobel x i y parametara za ulaznu sliku a
%   [sx,sy]=sobel_xy(a)
% V P, Februar 2000
% Calc. x and y parameters of the Sobel filter for the input image

% Definisemo Sobel filtere: Defining Sobel filter
tx=double([ -1 0 1; -2 0 2; -1 0 1]);
ty=double([ -1 -2 -1; 0 0 0; 1 2 1]);

% Definisemo velicinu rezultata: Defining the size of the result
if length(varargin)>0 & strcmp(varargin{1},'valid')
  velicina='valid';
else
  velicina='same';
end

% Filtriramo sobel filterima: Filtering Sobel
sx=filter2(tx,a,velicina)/8;
sy=filter2(ty,a,velicina)/8;

% Skidamo ivicne efekte ako se trazi ista velicina: Getting rid of the
% boundary effects if the same size or the result is required
if strcmp(velicina,'same');
  [v,s]=size(sx);
  sx(1,:)=zeros(1,s); sx(:,1)=zeros(v,1); sx(:,s)=zeros(v,1); sx(v,:)=zeros(1,s);
  sy(1,:)=zeros(1,s); sy(:,1)=zeros(v,1); sy(:,s)=zeros(v,1); sy(v,:)=zeros(1,s);
end