function [Y,L] = dht2(X,N,D,T,cod,varargin)
% DHT2 - 2-D discrete Hermite transform
%
% Synopsis:
% Y = DHT2(X,N,D,T)
% Y = DHT2(X,N,D,T,SHAPE)
% [Y,L] = DHT2(...,COD,P1,P2,...)
%
% Inputs:
% X - gray level image normalized in the interval [0 1]
% N - scale parameter (See DHTMTX)
% D - maximum order of the expansion (0 <= D <= 2*N)
% T - sampling distance (1 <= T <= N)
% COD - one of the following character that indicates an extra post-processing:
%     'r' - rotate coefficients (See RDHT)
%     'd' - directional-oriented coefficients (See DDHT)
%     'c' - classify rotated coefficients (See CRDHT)
%     's' - smooth coefficients (See SDHT)
%     'q' - use a quincunx lattice (See QDHT).
%     Use 'qr', 'qd', 'qc' or 'qs' to speed up the computation of
%     a post-processing over a quincunx lattice.
% SHAPE - (optional) use 'full', 'valid' or 'same' to control the size of
%     the coefficients as in CONV2. Use 'reflex' to reflex edges
%     before filtering the image and cropping to its original size.
% ... indicate extra parameters for post-processing
%
% Outputs:
% Y - transfrom coefficients ordered as described in DHTORD
% L - calss label for each lattice point assinged only if COD
%     contains one character 'c'.
%
% Used functions: DHTMTX, RDHT, DDHT, CRDHT, SDHT, QDHT
%
% See also IDHT2

% Support for multi-band images
K = size(X,3);
if K > 1,
    for k = 1:K,
        Y(:,:,:,k) = dht2(X(:,:,k),N,D,T);
    end
    return;
end

% Single-band images
shape = 'full';
if nargin < 5,
   cod = '';
elseif length(cod) > 2,
   shape = cod;
   if length(varargin),
      cod = varargin{1};
      varargin(1) = [];
   else
      cod = '';
   end   
end

if strcmp(shape,'reflex'),
   BC = ceil(N/2)+1:-1:2; BF = 1:floor(N/2);
   X = [X(BC,BC)     X(BC,:)     X(BC,end-BF)
        X(: ,BC)     X(: ,:)     X(: ,end-BF)
        X(end-BF,BC) X(end-BF,:) X(end-BF,end-BF)];
   shape = 'valid';
end

% Standard DHT

H = dhtmtx(N,D);  % Binomial functions
H(:,1);
i = 1;
t0 = clock;

for n = 0:min(D,2*N),
  for m = max(0,n-N):min(N,n),
      h1=H(:,m+1);
      h2=H(:,n-m+1);
      
    y = conv2(H(:,m+1),H(:,n-m+1),X,shape);
    Y(:,:,i) = y(1:T:end,1:T:end);
    i = i+1;
	end
end
if isempty(cod), return; end

%%%%%%%%%%%%%%%%%% Optional post-procesing %%%%%%%%%%%%%%%%%%

switch(lower(cod(1)))
case 'q'
   if length(cod) > 1,
      cod = cod(2);
   else
      cod = '';
   end
   % To avoid warnings
   if cod == 'c',
      [Y(:),L] = qdht(Y,N,D,cod,varargin{:});
   else
      Y(:) = qdht(Y,N,D,cod,varargin{:});
   end         
case 'r'
   Y(:) = rdht(Y,N,D,varargin{:});
case 'd'
   Y(:) = ddht(Y,N,D,varargin{:});
case 'c'
   [Y(:),L] = crdht(Y,N,D,varargin{:});
case 's'
   Y(:) = sdht(Y,N,D,varargin{:}); 
case 'g'
   offs = ceil(N/T);
   y = Y(1+offs:end-offs,1+offs:end-offs,1);
   Y(:,:,1) =  (1+tanh((Y(:,:,1)-mean2(y))/std2(y)))/2;   
case 'e'
   offs = ceil(N/T);
   Y(1+offs:end-offs,1+offs:end-offs,1) =  histeq(Y(1+offs:end-offs,1+offs:end-offs,1));   
end