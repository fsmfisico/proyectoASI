function [Z,Z0] = dhtshow(Y,N,D)
% DHTSHOW - show the 2-D coefficients
%
% Synopsis:
% DHTSHOW(Y,D) 
% DHTSHOW(Y,N,D)
%
% Inputs
% Y - matrix or cell array containing the 2-D transfrom coefficients 
% N - filters order (N+1 equals the filter length). This is the minimum scale
%     for multiscale decomposition 
% D - maximum order of the expansion
%
% See also DHT2, MDHT2, CLSSPLOT, ANGSHOW

  
if ~iscell(Y),
    if size(Y,4) == 3,
        for k = 1:3,
            Z(:,:,k) = min(max(dhtshow(Y(:,:,:,k),N,D),0),1);
        end
    else
        dims = [size(Y,1),size(Y,2)];
        Z = ones((min(N,D)+1)*dims);
        i = 1;
        for n = 0:min(2*N,D),
    %        c = (4/3)^n;
    c=(2)^n;
            for m = max(0,n-N):min(N,n),
                p = dims(1)*m+1:dims(1)*(m+1);  
                q = dims(2)*(n-m)+1:dims(2)*(n-m+1);
                Z(p,q) = c*Y(:,:,i)+0.5*(n>0);
                i = i+1;
            end
        end
    end
    if nargout > 0, return; end
    imshow(Z)
   clear Z
   return;
end

if nargin < 2 | ~iscell(Y), return; end

if nargin == 2,
   D = N;
end

% pyramid show
M = length(Y);
dims = size(Y{1})-2;
dims = dims(1:2);
twoplots = D > 2 & N == 2;
if twoplots,
  Z0 = ones(min(3,D+1)*dims);
  y = Y{1};
  y = y(2:end-1,2:end-1,:);
  N = 2;
  i = 1;
  for n = 1:D,
    c = (4/3)^n;
    for m = max(0,n-N):min(N,n),
      p = dims(1)*m+1:dims(1)*(m+1);  
      q = dims(2)*(n-m)+1:dims(2)*(n-m+1);
      Z0(p,q) = c*y(:,:,i);
      i = i+1;
    end
 end
 dims = dims/2;
end

Z = ones(min(7,D+1)*dims);

off = [0,0];
y = Y{M};
y0 = y(:,:,1);
Y{M} = y(:,:,2:end);
N = 6;

for k = 1+twoplots:M,   
  y = Y{k};
  y = y(k:k+dims(1)-1,k:k+dims(2)-1,:);
  i = 1;
  for n = 1:D,
    c = (4/3)^n;
    for m = max(0,n-N):min(N,n),    
      p = dims(1)*m+1+off(1):dims(1)*(m+1)+off(1);
      q = dims(2)*(n-m)+1+off(2):dims(2)*(n-m+1)+off(2);
        Z(p,q) = c*y(:,:,i);
%        Z(p,q) = mat2gray(y(:,:,i));
        i = i+1;
      end
  end
  dims = floor(dims/2);  
  off = off+(D+1)*dims;
end
dims = 2*dims;
y0 = y0(k:k+dims(1)-1,k:k+dims(2)-1)-0.5;
Z(end-dims(1)+1:end,end-dims(2)+1:end) = y0;

if nargout == 0,
   warning off
   if twoplots,
      figure; imshow(Z0+0.5);
   end   
   figure; imshow(Z+0.5)
   clear Z Z0
   warning on
end
