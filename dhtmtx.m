function [H,G] = dhtmtx(N,D,T)
% DHTMTX - discrete Hermite transform matrix
%
% Synopsis:
% H = DHTMTX(N) 
% H = DHTMTX(N,D)
% [H,G] = DHTMTX(N,D)
% [H,G] = DHTMTX(N,D,T)
% 
% Inputs:
% N - transform order. N+1 equals the filters length (N >= 0)
%     It has been shown in the literature that a limiting process 
%     turns the DHT into the continous Hermite transform,
%     which is based in gaussian derivatives of a given scale.
%     The scale parameter is tipically one half of the variance
%     of the Gaussian window, which is approximated by N/4.
% D - maximum order of the expanssion (0 <= D <= N)
% T - sampling distance (1 <= T <= N)
%
% Outputs:
% H - transform matrix so that H(:,k) is the filter function of order k-1
% G - interpolation matrix so that G(:,k) is the interpolating function of order k-1
%
% Used functions: CONV2, CONVMTX

% Last revision: 25/Feb/02 
% by J. L. Silvan-Cárdenas
% jlsilvan@centrogeo.org.mx

if nargin < 2,
  D = N;
else
  D = min(N,D);
end

if N == 0, H = 1; return; end

% Define the binomial masks
if D > 0,
  B = [1,1;1,-1];  
  %B = [sin(pi/6),cos(pi/6);cos(pi/6),-sin(pi/6)]; 
else
  %B = [1;1];
  B = [1;1];
end

H = B;
size(H,2);
for m = 2:N,
   % Compute the binomial filters by sucesive convolution
   if size(H,2) <= D,
      H = [conv2(B(:,1),1,H),conv2(B(:,2),1,H(:,end))];
   else
      H = conv2(B(:,1),1,H);
   end
end
C = sqrt(H(1:D+1,1))'/2^N; %Pa que sirve la raiz??
% C = sqrt(H(1:D+1,1))'/1;
% Identity=H*H
H(:) = H.*C(ones(N+1,1),:);

if nargout > 1 
   % Compute the interpolation functions
   if nargin < 3,       
      T = 1;
   end  
   W = convmtx(H(:,1),N+1);
  	W = sum(W(N+1-floor(N/T)*T:T:2*N+1,:))';
  	G = H(end:-1:1,:)./W(:,ones(1,D+1)); % pattern functions
end
