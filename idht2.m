function [X,x] = idht2(Y,xsiz,N,D,T,cod,varargin)
% IDHT2 - 2-D inverse discrete Hermite transform
%
% Synopsis: 
% X = IDHT2(Y,XSIZ,N,D,T)
% X = IDHT2(Y,XSIZ,N,D,T,'shape')
% [L,H] = IDHT2(Y,XSIZ,N,D,T)
%
% Inputs:
% Y - transform coefficients
% XSIZ - size of the original image
% N - transform order (See DHTMTX)
% D - maximum order of the expansion (0 <= D <= 2*N)
% T - sampling distance (1 <= T <= N)
% COD - one of the following character that indicates an extra pre-processing:
%     'r' - rotate coefficients (See RDHT)
%     'd' - directional-oriented coefficients (See DDHT)
%     'c' - classify rotated coefficients (See CRDHT)
%     's' - smooth coefficients (See SDHT)
%     'q' - use a quincunx lattice (See QDHT).
%     Use 'qr', 'qd', 'qc' or 'qs' to speed up the computation of
%     a post-processing over a quincunx lattice.
%     Use 'full', 'valid' or 'same' to control the size of
%     the coefficients as in CONV2.
% ... indicate extra parameters for pre-processing
%
% Outputs:
% X - resinthesized gray level image
% L - low-pass component of the signal
% H - high-pass component of the signal (X = L+H)
%
% Note: whenever you specify an extra-processing with a quincunx lattice
% and want to specify extra-parameters for the corresponding
% extra-processing function (i,e. RDHT) you must specify the parameter AC
% for the function QDHT first. 
% For example:
%
%       X = imread('rice.tif');
%       % rotate and classify coefficients over a quincunx lattice
%       [Y,L] = dht2(X,8,4,4,'qc'); % Here AC is 0 
%       % perform an inverse rotation and resynthesize the image
%       X2 = idht2(Y,size(X),8,4,4,'qr',0,'inv'); % pases AC=0 to QDHT which pases 'inv' to RDHT
%       subplot(1,2,1); imshow(X,[0,255]), title('original image')
%       subplot(1,2,2); imshow(X,[0,255]), title('resynthesized image')
%
% See also DHT2


K = size(Y,4);
if K > 1,
    for k = 1:K,
        X(:,:,k) = idht2(Y(:,:,:,k),xsiz(1:2),N,D,T);
    end
    return;
end

shape = 'valid';

if nargin < 6,
    cod = '';
elseif length(cod) > 2,
    shape = cod;
    cod = '';  
end

quincunx = 0;

%%%%%%%%%%%%%%%%%% Optional pre-processing %%%%%%%%%%%%%%%%%%
if ~isempty(cod),
    switch(lower(cod(1)))
        case 'q'
            quincunx = 1;
            T2 = 2*T;
            [H,G] = dhtmtx(N,D,T2);  % Hermite functions
            H(:) = flipud(H);
            i = rem((ceil((N+1)/2)-T*ceil(N/T2)+ceil(N/2):ceil((N+1)/2)+T2-T*ceil(N/T2)+ceil(N/2)-1)-1,N)+1;
            W = H(i,1)./G(i,1); W = W*W';
            W = repmat(W+fftshift(W),ceil(xsiz/T2)); % Weighting funtion    
            G(:) = H; 
            
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
            
            Y(isnan(Y)) = 0;
        case 'r'
            Y(:) = rdht(Y,N,D,varargin{:});
        case 'd'
            Y(:) = ddht(Y,N,D,varargin{:});
        case 'c'
            [Y(:),varargout{:}] = crdht(Y,N,D,varargin{:});
        case 's'
            Y(:) = sdht(Y,N,D,varargin{:}); 
        case 'g'
            offs = ceil(N/T);
            y = Y(1+offs:end-offs,1+offs:end-offs,1);
            Y(:,:,1) =  (1+tanh((Y(:,:,1)-mean2(y))/std2(y)))/2;   
    end
    if any(cod == 'e')
        Y(:) = edht(Y);
    end   
end   
%%%%%%%%%%%%%%%%%% Standard IDHT %%%%%%%%%%%%%%%%%%%%%%%

if ~quincunx,
    [H,G] = dhtmtx(N,D,T);  % Hermite functions     
end
switch(shape)
    case 'full', ysiz = xsiz+N; 
    case 'same', ysiz = xsiz; 
    case 'valid', ysiz = xsiz+N;
end
yi = zeros(ysiz);
yi(1:T:end,1:T:end) = Y(:,:,1);
X = conv2(G(:,1),G(:,1),yi,shape); % Low-pass component
G(:,1);%%%%%%%
x = 0; i = 1;
for n = 1:min(D,2*N),
    for m = max(0,n-N):min(N,n),
        i = i+1;    
        yi(1:T:end,1:T:end) = Y(:,:,i);
        g1=G(:,m+1);
        g2=G(:,n-m+1);
        x = x+conv2(G(:,m+1),G(:,n-m+1),yi,shape); % Hihg-pass component
    end
end

if quincunx,
    X(:) = (X+x)./W(1:xsiz(1),1:xsiz(2));     
elseif nargout == 1,
    X(:) = X+x;
end

