function z = ddht(y,N,D,flag,theta)
% DDHT - direct the 2-D coefficients
%
% Synopsis:
% Yd = DDHT(Yc,N,D)
% Yd = DDHT(Yc,N,D,DIR)
% Yc = DDHT(Yd,N,D,DIR,THETA)
%
% Inputs:
% Yc - cartesian coefficients 
% Yd - directional-oriented coefficients 
% N  - transform order
% D  - maximum order for the expansion
% DIR - either 'fwd' or 'inv' to specify the rotation 
%       direction (default 'fwd')
% THETA - Rotation angle. By default THETA is set to 
%      the gradient angle and Yd(:,:,2) = THETA/(2*PI)
%
% Outputs:
% Yc - cartesian coefficients 
% Yd - directional-oriented coefficients 
%
% See also: RDHT, SDHT, CRDHT

if nargin < 4,
   flag = 'fwd';
end
ysiz = size(y);
z = y;

if length(ysiz) < 3,
   return;
end

switch(flag)
case 'inv',
   
   if nargin < 5
   	theta = (2*pi)*y(:,:,3);
     	y(:,:,3) = 0;
     	c = cos(theta); s = sin(theta);
     	z(:,:,2) = c.*y(:,:,2);
     	z(:,:,3) = s.*y(:,:,2);
   else
      if length(theta) == 1
         theta = theta(ones(ysiz(1:2)));
      end
      c = cos(theta); s = sin(theta);
      z(:,:,2) = c.*y(:,:,2)-s.*y(:,:,3);            
      z(:,:,3) = c.*y(:,:,3)+s.*y(:,:,2);      
   end
     i = 2:3;  
  % steering of coeffitients above the main diagonal  
  for n = 2:min(D,N),
    i = i(end)+1:i(end)+n+1;
    npi = shiftdim((0:n)*(pi/(n+1)),-1); 
    thetaj = theta(:,:,ones(1,1,n+1))+npi(ones(1,ysiz(1)),ones(1,ysiz(2)),:);
    c = cos(thetaj);  s = sin(thetaj);
    for m = 0:n,
      z(:,:,i(m+1)) = sum(y(:,:,i).*(c.^(n-m)).*(s.^m),3);
    end
    an = sum(sin((0:n)*(pi/(n+1))).^(2*n));
    C = shiftdim(sqrt(poly(-ones(1,n)))/an,-1);
    z(:,:,i) = C(ones(1,ysiz(1)),ones(1,ysiz(1)),:).*z(:,:,i);
  end
    
  % steering of coeffitients below the main diagonal  
  for n = 1:min(D-N,N-1),
    nn = N-n;
    i = i(end)+1:i(end)+nn+1;
    thetaj = repmat(theta,[1,1,nn+1])+repmat(shiftdim(pi*(0:nn)/(nn+1),-1),[ysiz(1:2),1]);
    an = sum(sin(thetaj).^(2*nn),3);
    c = cos(thetaj);  s = sin(thetaj);
    for m = 0:nn,
      C = sqrt(prod(nn-m+1:nn)/prod(1:m));
      z(:,:,i(m+1)) = sum(y(:,:,i).*C.*(c.^(nn-m)).*(s.^m),3)./an;
    end
  end

case 'fwd'
   
   if nargin < 5
      theta = atan2(y(:,:,3),y(:,:,2)); % gradient angle
      c = cos(theta); s = sin(theta);      
      z(:,:,2) = c.*y(:,:,2)+s.*y(:,:,3);      
      z(:,:,3) = theta/(2*pi);
   else
      if length(theta) == 1
         theta = theta(ones(ysiz(1:2)));
      end      
      c = cos(theta); s = sin(theta);
      z(:,:,2) = c.*y(:,:,2)+s.*y(:,:,3);            
      z(:,:,3) = c.*y(:,:,3)-s.*y(:,:,2);
   end    
   i = 2:3;
  	% steering of coeffitients above the main diagonal
  	for n = 2:min(D,N),
   	i = i(end)+1:i(end)+n+1;
    	dtheta = pi/(n+1);  % angle increment
    	m = repmat(shiftdim(0:n,-1),[ysiz(1:2),1]);
    	C = repmat(shiftdim(sqrt(poly(-ones(1,n))),-1),[ysiz(1:2),1]);
    	for j = 0:n,
       	c = repmat(cos(theta+j*dtheta),[1,1,n+1]);
       	s = repmat(sin(theta+j*dtheta),[1,1,n+1]);
       	z(:,:,i(j+1)) = sum(y(:,:,i).*C.*(c.^(n-m)).*(s.^m),3);
    	end
  	end
  
  	% steering of coeffitients below the main diagonal
  	for n = 1:min(D-N,N-1),    
    	nn = N-n;
    	i = i(end)+1:i(end)+nn+1; 
    	dtheta = pi/(nn+1);  % angle increment
    	m = repmat(shiftdim(0:nn,-1),[ysiz(1:2),1]);
    	W = repmat(shiftdim(sqrt(poly(-ones(1,nn))),-1),[ysiz(1:2),1]);
    	for j = 0:nn,
       	c = repmat(cos(theta+j*dtheta),[1,1,nn+1]);
       	s = repmat(sin(theta+j*dtheta),[1,1,nn+1]);
       	z(:,:,i(j+1)) = sum(y(:,:,i).*W.*(c.^(nn-m)).*(s.^m),3);
    	end
 	end

end

