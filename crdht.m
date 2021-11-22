function [Y,F] = crdht(Y,N,D,T,fact)
% CRDHT - rotatate and classify the 2-D coefficients
%
% Synopsis:
% [Yc,L] = CRDHT(Y,N,D)
% [Yc,L] = CRDHT(Y,N,D,THR)
%
% Inputs:
% Y - transfrom coefficients
% N - transfrom order
% D - maximum order of the expansion (0 <= D <= 2*N)
% THRE - vector used in the classification (default is [0.2, 0.8])
%
% Outputs:
% Yc - classified transfrom coefficients. Depending on the class they belong to,
%      irrelevant coefficients are set to zero.
% L - Class labels. 0 for 0-D patterns, 1 for 1-D patterns, 2 for 2-D patterns.
%
% Note: This function is based on the method described in 
% in Chap. 3 "Clasificación de estructuras en el dominio de la THDR"
% of the master tesis by J. Luis Silván-Cárdenas 
% from Div. de Estudios de Posgrado, Fac. de Ing., UNAM, 2002, 
% entitled "Compresión de imágenes basada en modelos 
% Gaussianos de percepción visual"
%
% Used functions: RDHT
%
% See also: DDHT, SDHT, RDHT

if nargin < 4, fact = [0.1, 0.9]; end

ysiz = size(Y);
Y(:) = rdht(Y,N,D); % Perform the rotation

L = Y(:,:,1);

i = 2; k = 1;
dE = 0;
for n = 1:min(D,2*N),
  dE = dE+Y(:,:,i).^2;
  for m = max(0,n-N):min(N,n),
    if m > 0, idx(k) = i; k = k+1; end
    i = i+1;
  end
end
idx = idx(2:end); % drop the index for the angle

E = sum(Y(:,:,[2,4:end]).^2,3);
dE(:) = sqrt(E-dE);
E(:) = sqrt(E);

a = 0.7; b = 0.1; c = 0.35; d = 0.8;

% Separation of class 0-D
thr0 = fact(1)*(b+abs(L.^a-c^a).^(1/a)./(L.^a+c^a).^(1/a));
F = E > thr0;
Y(:,:,2:end) = Y(:,:,2:end).*F(:,:,ones(1,1,ysiz(3)-1));

% Separation of class 1-D
thr1 = fact(2)*max(thr0,E.^d.*thr0.^(1-d));
F1 = (F & dE > thr1);
Y(:,:,idx) = Y(:,:,idx).*F1(:,:,ones(1,1,length(idx)));
F(:) = F+F1;
