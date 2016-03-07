function varargout = PolyDeriv(Y, X, NumPoints, Order)
% [dY_1, dY_2, ...] = PolyDeriv(Y, DeltaX, NumPoints, Order)
%    ... OR ...
% [dY_1, dY_2, ...] = PolyDeriv(Y, X, NumPoints, Order)
% Takes numerical derivatives by fitting polynomials to Y.
% If X has constant spacing, using the first form is
% computationally more efficient.
%  INPUT:
%   -Y:  waveform to be differentiated
%   -DeltaX: sampling period
%   -X: list of sample locations/times
%   OPTIONAL:
%    -NumPoints: number of points to use in fitting polynomials
%      (default is 5)
%    -Order: order of polynomial to use in fitting. (default is 
%      1 + maximum requested derivative if NumPoints is specified,
%      otherwise it is 4)
%  OUTPUT:
%   -dY_1:  The first derivative of Y.  It is sampled at the same X
%     coordinates as Y.
%   ... each subsequent output argument is a higher-order
%   derivative.  You can request at most Order derivatives.


if(nargin < 2)
  error('Incorrect number of input arguments')
elseif(nargin < 3)
  NumPoints = 5;
  Order = 4;
elseif(nargin < 4)
  Order = nargout + 1;
  if(NumPoints ~= round(NumPoints))
    error('NumPoints must be an integer')
  end
else
  if(NumPoints ~= round(NumPoints))
    error('NumPoints must be an integer')
  elseif(Order ~= round(Order))
    error('Order must be an integer')
  end
end
if(nargout > Order)
  error(['Too many output arguments.  To take higher derivatives,' ...
	 ' increase Order'])
end
if(mod(NumPoints, 2) == 0)
  error('NumPoints must be odd')
end

if(size(Y,2) > size(Y,1))
  Y = Y';
  Flip = true;
else
  Flip = false;
end

if(length(X) == 1)
  DeltaX = X;
  Coefs = GetDerivCoefs(NumPoints, Order, DeltaX);
  NumLeadPts = (NumPoints - 1) / 2;
  if(nargout == 1)
    YPrime = zeros(size(Y));
    TempYFront = Y(1:NumPoints,:);
    TempYBack = Y((end-NumPoints+1):end,:);
    
    for n = 1:NumLeadPts
      C1 = Coefs{n}(1,:);
      C2 = Coefs{NumPoints + 1 - n}(1,:);
      YPrime(n,:) = C1 * TempYFront;
      YPrime(end-n+1,:) = C2 * TempYBack;
    end
    
    C = Coefs{NumLeadPts + 1}(1,:);
    Temp = zeros(size(Y,1)-NumPoints+1, size(Y,2));
    for n = 1:NumPoints
      Temp = Temp + C(n) * Y(n:(end-NumPoints+n),:);
    end
    YPrime(NumLeadPts+1:(end-NumLeadPts),:) = Temp;
    
    varargout = {YPrime};
  else
    varargout = cell(nargout, 1);
    Temp = cell(nargout, 1);
    for DOrder = 1:nargout
      varargout{DOrder} = zeros(size(Y));
      Temp{DOrder} = zeros(size(Y,1)-NumPoints+1, size(Y,2));
    end
    TempYFront = Y(1:NumPoints,:);
    TempYBack = Y((end-NumPoints+1):end,:);
    
    for n = 1:NumLeadPts
      C1 = Coefs{n};
      C2 = Coefs{NumPoints + 1 - n};
      for DOrder = 1:nargout
	varargout{DOrder}(n,:) = C1(DOrder,:) * TempYFront;
	varargout{DOrder}(end-n+1,:) = C2(DOrder,:) * TempYBack;
      end
    end
    
    C = Coefs{NumLeadPts + 1};
    for n = 1:NumPoints
      TempY = Y(n:(end-NumPoints+n),:);
      for DOrder = 1:nargout
	Temp{DOrder} = Temp{DOrder} + C(DOrder,n) * TempY;
      end
    end
    for DOrder = 1:nargout
    varargout{DOrder}(NumLeadPts+1:(end-NumLeadPts),:) = Temp{DOrder};
    end
  end
else  % X is a list of values
  if(size(X,1) == 1)
    X = X';
  end
  if(length(X) ~= length(Y))
    error('X and Y must have the same number of sample points.')
  end
  varargout = cell(nargout, 1);
  for DOrder = 1:nargout
    varargout{DOrder} = zeros(size(Y));
  end
  TempYFront = Y(1:NumPoints,:);
  TempYBack = Y((end-NumPoints+1):end,:);
  TempXFront = X(1:NumPoints);
  TempXBack = X((end-NumPoints+1):end);
  
  NumLeadPts = (NumPoints - 1) / 2;
  for n = 1:NumLeadPts
    C1 = GetInstDerivCoefs(NumPoints, Order, TempXFront - X(n));
    C2 = GetInstDerivCoefs(NumPoints, Order, TempXBack - X(end-n+1));
    for DOrder = 1:nargout
      varargout{DOrder}(n,:) = C1(DOrder,:) * TempYFront;
      varargout{DOrder}(end-n+1,:) = C2(DOrder,:) * TempYBack;
    end
  end
  
  nStart = NumLeadPts + 1;
  nStop = length(X) - NumLeadPts;
  for n = nStart:nStop
    TempX = X((n-NumLeadPts):(n+NumLeadPts));
    TempY = Y((n-NumLeadPts):(n+NumLeadPts),:);
    C = GetInstDerivCoefs(NumPoints, Order, TempX - X(n));
    for DOrder = 1:nargout
      varargout{DOrder}(n,:) = C(DOrder,:) * TempY;
    end
  end
end

if(Flip)
  for n=1:nargout
    varargout{n} = varargout{n}';
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Coefs = GetDerivCoefs(NumPoints, Order, DeltaX)

J = zeros(Order + 1, NumPoints);
for n=1:NumPoints
  Z = (1-n):(NumPoints-n);
  for m=0:Order
    J(m+1,:) = Z.^m;
  end
  C = pinv(J');
  C = C(2:end,:);
  Fact = 1;
  for m=1:Order
    Fact = Fact * m / DeltaX;
    C(m,:) = C(m,:) * Fact;
  end
  Coefs{n} = C;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = GetInstDerivCoefs(NumPoints, Order, X)

JTrans = zeros(NumPoints, Order + 1);

for m=0:Order
  JTrans(:,m+1) = X.^m;
end
C = pinv(JTrans);
C = C(2:end,:);
Fact = 1;
for m=2:Order
  Fact = Fact * m;
  C(m,:) = C(m,:) * Fact;
end

return