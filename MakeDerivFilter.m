function dFilt = MakeDerivFilter(order, dx, fPass, fStop)
% dFilt = MakeDerivFilter(order, dx, fPass, fStop)
% Constructs a filter that takes noise-filtered derivatives of the
% requested order.
%
% Designed to  be used in conjunction with FilterSignal.m
%
% Note that the quality of the filtered signal will be lower near
% the beginning and end.
%  INPUTS:
%   -order: order of the derivative
%   -dx: sample spacing
%   -fPass: low-pass frequency (this frequency and below should
%           pass through the filter)
%   -fStop: high-stop frequency (this frequency and above should be
%           stopped by the filter)
%  OUTPUTS:
%   -dFilt: length-L cell array, When passed to FilterSignal, it
%           will produce an array of derivatives sampled at the
%           same times as the original signal.

if nargin ~= 4
  help MakeDerivFilter
  error('Invalid number of inputs')
end
if nargout ~= 1
  help MakeDerivFilter
  error('Invalid number of outputs')
end

wPass = 2*pi*dx*fPass;
wStop = 2*pi*dx*fStop;

halfLen = round(1.0 / (dx * fPass));
filterLen = 2 * halfLen + 1;
dFilt = cell(1, filterLen);

for delay=-halfLen:halfLen
  n = delay + halfLen + 1;
  dFilt{n} = constructFilter(order, dx, filterLen, delay, wPass, wStop);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dFilt = constructFilter(order, dx, filterLen, delay, ...
				 wPass, wStop)
useShift = false;
if ~useShift
  minLen = max(round(0.5 * pi / wPass), order + mod(order, 2) + 5);
  filterLen = filterLen - 2 * abs(delay);
  if filterLen < minLen
    delay = sign(delay) * round((minLen - filterLen) / 2);
    filterLen = minLen;
    % Can't (or don't know how) use this technique to make a
    % reliable filter when the filterLen is this short.  So make a
    % polynomial (Savitzky-Golay) filter instead.  It will at least
    % produce sensible results, even if they aren't optimal.
    polyOrder = max(order + 1, 2);
    dFilt = getPolyFilter(polyOrder, order, dx, filterLen, delay);
    return
  else
    delay  = 0;
  end
end

%The coefficients of this filter c_k are solutions to the equations
% For frequency w in the pass band:
%  sum_k=-m:m c_k e^(i w k) = e^(i w delay) (i w)^order
% For frequency w in the stop band:
%  sum_k=-m:m c_k e^(i w k) = 0

%To make the equations overdetermined, use roughly twice as many
% equations as needed (use real and imaginary parts, and treat as
% separate equations)
%Use half for the pass band, sampling evenly from w in [0 wPass]
% and half in the stop band, sampling evenly from w in [wStop pi]
halfLen = (filterLen - 1) / 2;
numEq = (halfLen + 1) * 4;
numPass = halfLen + 1;
numStop = numPass;

mat = zeros(numEq, filterLen);
vec = zeros(numEq, 1);

orderEven = (mod(order, 2) == 0);

fourierInd = -halfLen:halfLen;
wVec = linspace(0, wPass, numPass);
rowReal = -1;
for m = 1:numPass;
  rowReal = rowReal + 2;
  rowImag = rowReal + 1;
  
  w = wVec(m);
  mat(rowReal,:) = cos(w*fourierInd);
  mat(rowImag,:) = sin(w*fourierInd);
  
  if orderEven
    w_order = w^order * (-1)^(order/2);
    vec(rowReal) = w_order * cos(w * delay);
    vec(rowImag) = w_order * sin(w * delay);
  else
    w_order = w^order * (-1)^((order - 1)/2);
    vec(rowImag) = w_order * cos(w * delay);
    vec(rowReal) = -w_order * sin(w * delay);
  end
end
wVec = linspace(wStop, pi, numStop);
for m = 1:numStop
  rowReal = rowReal + 2;
  rowImag = rowReal + 1;
  
  w = wVec(m);
  mat(rowReal,:) = cos(w*fourierInd);
  mat(rowImag,:) = sin(w*fourierInd);
end

%calculate the dx-independent coeffients of the filter:
dFilt = pinv(mat) * vec;

%reverse for convolution, multiply by appropriate power of dx:
dFilt = dFilt(filterLen:-1:1) * (dx^-order);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dFilt = getPolyFilter(polyOrder, order, dx, filterLen, ...
					  delay)
% At edges, fall back to using Savitzky-Golay filter, which is
% reliable even if it's difficult to precisely tune its frequency
% characteristics.
J = zeros(filterLen, polyOrder + 1);
halfLen = (filterLen - 1) / 2;
n = (-halfLen-delay):(halfLen-delay)';

for m = 0:polyOrder
  J(:, m+1) = n.^m;
end
JInv = pinv(J);

%calculate the dx-independent coeffients of the filter:
dFilt = factorial(order) * JInv(order+1,:);

%reverse for convolution, multiply by appropriate power of dx:
dFilt = dFilt(filterLen:-1:1) * (dx^-order);
return
