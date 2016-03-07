function [ySpectrum, varargout] = Spectrum(y, dT, varargin)
% [ySpectrum, f, yCorr] = Spectrum(y, dt, method)

% get options
defaultOptions = { ...
  'removeTrend', true, ...
  'plotSubject', false ...
};
options = GetOptions(defaultOptions, varargin);

if size(y,1) > 1
  y = y';
end

if options.removeTrend
  % remove linear trend from y, unless requested not to
  y = removeTrend(y, dT);
end

% compute spectrum
yCorr = xcorr(y, 'unbiased');
numY = length(y);
halfInd = ceil(numY / 2);
ySpectrum = fft(yCorr(halfInd:(halfInd + numY - 1)));

% determine if a plot is needed
needPlot = ischar(options.plotSubject) || options.plotSubject;

if needPlot || nargout > 1
  % need to compute the frequencies
  numSpec = length(ySpectrum);
  halfInd = ceil(numSpec/2);
  f = (2*pi / numSpec) * (0:(numSpec-1));
  f(halfInd:end) = f(halfInd:end) - 2*pi;
  if nargin > 1 && ~isempty(dT)
    f = f / dT;
  end
  if needPlot
    titleStr = 'Spectrum';
    if ischar(options.plotSubject)
      titleStr = [options.plotSubject, ' - ', titleStr];
    end
    h = NamedFigure(titleStr);
    set(h, 'WindowStyle', 'docked')
    clf
    hold on
    plot(f(1:(halfInd-1)), sqrt(abs(ySpectrum(1:(halfInd-1)))))
    plot(f(halfInd:end), sqrt(abs(ySpectrum(halfInd:end))))
    hold off
  end
  
  if nargout > 1
    varargout = {f};
  end
else
  varargout = {};
end

if nargout > 2
  varargout = [varargout, yCorr];
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = removeTrend(y, dt)
% remove linear trend from y
numY = length(y);
tEnd = dt * (numY - 1);
t = 0:dt:tEnd;
dYAvg = (y(end) - y(1)) / tEnd;
tMean = tEnd / 2;
yAvg = mean(y) + dYAvg * (t - tMean);
y = y - yAvg;
return