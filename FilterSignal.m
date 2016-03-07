function yFiltered = FilterSignal(y, filt)
% yFiltered = FilterSignal(y, filt)
% Filters a signal.  This function was designed to be used in
% conjuction with MakeDerivFilter.m
%  INPUTS:
%   -y: length-N vector, signal to be filtered
%   -filt: length-L cell array.  L must be less than N
%          each element is a filter vector that will be convolved
%          with y
%  OUTPUTS:
%   -yFiltered: length-N vector, filtered signal
if nargin ~= 2
  help FilterSignal
  error('Invalid number of inputs')
end
if nargout ~= 1
  help FilterSignal
  error('Invalid number of outputs')
end

filterLen = length(filt);
nHalf = (filterLen - 1) / 2;
yFiltered = zeros(size(y));

for n = 1:nHalf
  length_n = length(filt{n});
  yFront = y(1:length_n);
  yBack = y((end-length_n+1):end);
  yFiltered(n) = conv(yFront, filt{n}, 'valid');
  yFiltered(end-n+1) = conv(yBack, filt{end-n+1}, 'valid');
end

yFiltered((nHalf+1):(end-nHalf)) = conv(y, filt{nHalf + 1}, 'valid');
return