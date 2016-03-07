function [dy, varargout] = DerivFilter(y, dx, fPass, fStop)
% [dy, varargout] = DerivFilter(y, dx, fPass, fStop)
% Filter 
tPass = 1.0 / fPass;
tStop = 1.0 / fStop;
lenPass = tPass / dx;
lenStop = tStop / dx;
wPass = 2*pi / lenPass;
wStop = 2*pi / lenStop;

filterLen = 2 * round(lenPass) + 1;

dFilt = constructFilterNoDelay(1, dx, filterLen, wPass, wStop);
dy = applyFiltNoDelay(y, dFilt);

varargout = cell(1, nargout - 2);
for order = 2:nargout
  dFilt = constructFilterNoDelay(order, dx, filterLen, wPass, wStop);
  dyOrder = applyFiltNoDelay(y, dFilt);
  varargout{order-1} = dyOrder;
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dy, varargout] = DerivFilterOld(y, dx, fPass, fStop)

tPass = 1.0 / fPass;
tStop = 1.0 / fStop;
lenPass = tPass / dx;
lenStop = tStop / dx;
wPass = 2*pi / lenPass;
wStop = 2*pi / lenStop;

filterLen = 4 * round(lenPass) + 1;
halfLen = (filterLen - 1) / 2;
dFilt = cell(1, filterLen);

for delay=-halfLen:halfLen
  n = delay + halfLen + 1;
  dFilt{n} = constructFilter(1, dx, filterLen, delay, wPass, wStop);
end
dy = applyFilt(y, dFilt);

varargout = cell(1, nargout - 2);
for order = 2:nargout
  for delay=-halfLen:halfLen
    n = delay + halfLen + 1;
    dFilt{n} = constructFilter(order, dx, filterLen, delay, wPass, wStop);
  end
  d_order_y = applyFilt(y, dFilt);
  varargout{order-1} = d_order_y;
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dFilt = constructFilter(order, dx, filterLen, delay, wPass, wStop)
useShift = true;
if ~useShift
  minLen = max(round(0.5 * pi / wPass), order + mod(order, 2) + 5);
  filterLen = filterLen - 2 * abs(delay);
  if filterLen < minLen
    delay = sign(delay) * round((minLen - filterLen) / 2);
    filterLen = minLen;
    % Can't (or don't know how) use this technique to make a
    % reliable filter when the filterLen is this short.  So make a
    % polynomial filter instead.  It will at least produce sensible
    % results, even if they aren't optimal.
    polyOrder = max(order + 1, 2);
    dFilt = getPolyFilter(polyOrder, order, dx, filterLen, delay);
    return
  else
    delay  = 0;
  end
end

% create list of frequencies where filter behavior is specified (pass or
% stop)
halfLen = (filterLen - 1) / 2;
numPass = halfLen + 1;
numStop = halfLen;
wPassVec = linspace(0, wPass, numPass)';
wStopVec = linspace(wStop, pi, numStop)';

% create mat and vec to hold linear equation for filter coefficients
%mat = zeros(numRealEq, filterLen);
%vec = zeros(numRealEq, 1);

% fill out mat and vec according to linear equation
%   for each w in wPassVec
%     sum n=-hL:hL {e^(i n w)} = (iw)^order e^(i w delay)
%   for each w in wStopVec
%     sum n=-hL:hL {e^(i n w)} = 0
% note that equations are complex values, so break them up into pairs of
%   real-valued equations. Thus, divide matrix rows into 4 separate blocks:
%     real pass, imaginary pass, real stop, imaginary stop
n = -halfLen:halfLen;
mat = [cos(wPassVec * n); sin(wPassVec * n); ...
       cos(wStopVec * n); sin(wStopVec * n)];

orderEven = (mod(order, 2) == 0);
if orderEven
  w_order = wPassVec.^order * (-1)^(order/2);
  vec = [w_order .* cos(wPassVec * delay); ...
         w_order .* sin(wPassVec * delay) ; ...
         zeros(2 * numStop, 1)];
else
  w_order = wPassVec.^order * (-1)^((order - 1)/2);
  vec = [-w_order .* sin(wPassVec * delay); ...
         w_order .* cos(wPassVec * delay) ; ...
         zeros(2 * numStop, 1)];
end

%calculate the dx-independent coeffients of the filter:
dFilt = linsolve(mat, vec);

%reverse for convolution, multiply by appropriate power of dx:
dFilt = dFilt(filterLen:-1:1) * (dx^-order);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dFilt = constructFilterNoDelay(order, dx, filterLen, wPass, wStop)
% compute filter convolution coefficients

% create list of frequencies where filter behavior is specified (pass or
% stop)
halfLen = (filterLen - 1) / 2;
numPass = halfLen + 1;
numStop = halfLen;
wPassVec = linspace(0, wPass, numPass)';
wStopVec = linspace(wStop, pi, numStop)';

% create mat and vec to hold linear equation for filter coefficients
%mat = zeros(numRealEq, filterLen);
%vec = zeros(numRealEq, 1);

% fill out mat and vec according to linear equation
%   for each w in wPassVec
%     sum n=-hL:hL {e^(i n w)} = (iw)^order e^(i w delay)
%   for each w in wStopVec
%     sum n=-hL:hL {e^(i n w)} = 0
% note that equations are complex values, so break them up into pairs of
%   real-valued equations. Thus, divide matrix rows into 4 separate blocks:
%     real pass, imaginary pass, real stop, imaginary stop
n = -halfLen:halfLen;
mat = [cos(wPassVec * n); sin(wPassVec * n); ...
       cos(wStopVec * n); sin(wStopVec * n)];

orderEven = (mod(order, 2) == 0);
if orderEven
  w_order = wPassVec.^order * (-1)^(order/2);
  vec = [w_order; zeros(numPass + 2 * numStop, 1)];
else
  w_order = wPassVec.^order * (-1)^((order - 1)/2);
  vec = [zeros(numPass, 1); w_order; zeros(2 * numStop, 1)];
end

%calculate the dx-independent coeffients of the filter:
dFilt = linsolve(mat, vec);

%reverse for convolution, multiply by appropriate power of dx:
dFilt = dFilt(filterLen:-1:1) * (dx^-order);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dFilt = getPolyFilter(polyOrder, order, dx, filterLen, delay)
J = zeros(polyOrder + 1, filterLen);
halfLen = (filterLen - 1) / 2;
n = (-halfLen-delay):(halfLen-delay);

for m = 0:polyOrder
  J(m+1,:) = n.^m;
end
JInv = pinv(J');

%calculate the dx-independent coeffients of the filter:
dFilt = factorial(order) * JInv(order+1,:);

%reverse for convolution, multiply by appropriate power of dx:
dFilt = dFilt(filterLen:-1:1) * (dx^-order);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffy = applyFilt(y, dFilt)
% apply filter (calculating begining and ends of signal by using special
%  filters that read out shifted values)
% NOTE: this is VERY EXPENSIVE computationally, for crummy results
fLen = length(dFilt);
nHalf = (fLen - 1) / 2;
diffy = zeros(size(y));

for n = 1:nHalf
  fLen = length(dFilt{n});
  yFront = y(1:fLen);
  yBack = y((end-fLen+1):end);
  diffy(n) = conv(yFront, dFilt{n}, 'valid');
  diffy(end-n+1) = conv(yBack, dFilt{end-n+1}, 'valid');
end
diffy((nHalf+1):(end-nHalf)) = conv(y, dFilt{nHalf + 1}, 'valid');
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffy = applyFiltNoDelay(y, dFilt)
% apply filter, set the ends to be constant
fLen = length(dFilt);
nHalf = (fLen - 1) / 2;
diffy = zeros(size(y));

diffy((nHalf+1):(end-nHalf)) = conv(y, dFilt, 'valid');
diffy(1:nHalf) = diffy(nHalf + 1);
diffy(end-nHalf+1:end) = diffy(end-nHalf);
return