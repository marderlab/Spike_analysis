function [T, V, varargout] = CleanDcc(t, v, varargin)
% [T, V, Arr1, Arr2, ..., DccInfo] = CleanDcc(t, v, arr1, arr2, ...)
% Resamples a voltage trace to a courser sample, in order to reduce
% noise introduced by DCC.
% DccFreq is calculated by looking for the dominant fast frequency
% of abs(dv/dt)  (because voltage shouldn't change much during the
% DCC period).
%  INPUTS:
%   -t: Array of times
%   -v: Num_t x NumTraces array of voltages
%   -arr1: (OPTIONAL) Num_t x NumTraces array of voltages/currents.
%  OUTPUTS:
%   -T: Array of times
%   -V: Num_T x NumTraces array of voltages
%   -arr1: (OPTIONAL) Num_T x NumTraces array of currents, only
%       passed back if arr1 is passed in as non-empty array.
%   -DccInfo: (OPTIONAL) A structure containing DCC information
%      .DccFreq:  The frequency of DCC sampling (in kHz)
%      .DccPower: The spectrum power at the dccFreq
%      .f:         An array of spectrum frequencies
%      .Power:     An Array of spectrum powers

debugClean = false;

if nargin < 2
  error('Incorrect number of input arguments')
end
deltaNumArgs = nargout - nargin;
if deltaNumArgs < 0 || deltaNumArgs > 1
  help CleanDCC
  error('Invalid combination of in/out arguments.')
end

% First calculate the derivative and fast-frequency correlogram
numAutoCorr = find(t - t(1) > 50, 1);  %Specify length, 50 ms worth
[autoCorr, dV] = getAutoCorr(v, numAutoCorr);

% Next find the first maximum of autocorrelogram (First guess for DCC
%  period)
[dccInd, dccCorr] = getAutoCorrMax(autoCorr, debugClean, t);

if isnan(dccInd)  %Couldn't find dominant fast frequency, so no DCC
  dccFreq = Inf;
else
  % refine the estimation of the dcc power using spectral methods
  [dccFreq, dccPower, f, power] = ...
    refineDccFreq(t, dV, dccInd, dccCorr, debugClean);
end

varargout = cell(1, nargout - 2);
if isfinite(dccFreq)
  %Interpolate T, V, and if necessary, I
  T = t(1):(1/dccFreq):t(end);
  V = interp1(t, v, T);
  for n = 1:length(varargin)
    varargout{n} = interp1(t, varargin{n}, T);
  end
else
  T = t;
  V = v;
  for n = 1:length(varargin)
    varargout{n} = varargin{n};
  end  
end

if deltaNumArgs == 1
  dccInfo.dccFreq = dccFreq;
  dccInfo.dccPower = dccPower;
  dccInfo.f = f;
  dccInfo.power = power;
  varargout{nargout - 2} = dccInfo;
end

if debugClean
  h = NamedFigure('Voltage Trace'); %#ok<UNRCH>
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(t, v, 'r.', 'MarkerSize', 6)
  hold on
  plot(T, V, 'b.', 'MarkerSize', 6)
  hold off
  %keyboard
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [autoCorr, dV] = getAutoCorr(v, numAutoCorr)
% calculate the derivative and fast-frequency correlogram

% set dV to be the square of changes in v.
dV = diff(v).^2;
% Note that v may have many traces (v of form numSamples x numTraces)
if size(dV, 2) > 1
  % get the max dV for each trace
  maxDvByTrace = max(dV);
  % get the trace with the largest max dV
  [~, traceInd] = max(maxDvByTrace);
else
  traceInd = 1;
end
% set dV as the zScore of the trace with the largest max dV
dV = zscore(dV(:,traceInd));

% if there are too many samples, look at a subset of the trace
maxLen = 2^20 - 1;
if length(dV) > maxLen
  dV = dV((end-maxLen+1):end);
end

% calculate the autocorrelogram of this trace, restricted to the number of
%  requested elements in the autocorrelogram
autoCorr = xcorr(dV, numAutoCorr - 1, 'unbiased');
autoCorr = autoCorr(numAutoCorr:end);

% scale the autocorrelogram to start at 1.0
autoCorr = autoCorr / autoCorr(1);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dccInd, dccCorr] = getAutoCorrMax(autoCorr, debugClean, t)
% Find the first maximum of autocorrelogram (First guess for DCC period)

fSample = 1.0 / (t(2) - t(1));
[corrPower, corrFreq] = pmtm(autoCorr, 7/2, [], fSample);
ind1 = find(corrFreq >= 0.25, 1);
ind2 = find(corrFreq > 10, 1) - 1;
if isempty(ind2)
  ind2 = length(corrFreq);
end

[~, maxInd] = max(corrPower(ind1:ind2));
maxInd = maxInd + ind1 - 1;
dccInd = 1 + round(fSample / corrFreq(maxInd));
if dccInd == 1
  dccInd = 2;
elseif dccInd == length(autoCorr)
  dccInd = dccInd - 1;
end
if autoCorr(dccInd + 1) > autoCorr(dccInd)
  dccInd = dccInd + 1;
elseif autoCorr(dccInd - 1) > autoCorr(dccInd)
  dccInd = dccInd - 1;
end
dccCorr = autoCorr(dccInd);

if debugClean
  h = NamedFigure('CleanDCC Autocorrelation');
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(t(1:length(autoCorr))-t(1), autoCorr);
  hold on
  plot(t(dccInd)-t(1), dccCorr, 'ro', 'MarkerFaceColor', 'r')
  hold off
  xlabel('Time (ms)', 'FontSize', 18);
  ylabel('Autocorrelation', 'FontSize', 18);
  title('CleanDCC Autocorrelation', 'FontSize', 18);
  
  h = NamedFigure('CleanDCC Autocorrelation Power');
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(corrFreq, corrPower)
  xlabel('Frequency (kHz)', 'FontSize', 18);
  ylabel('Power', 'FontSize', 18);
  title('CleanDCC Autocorrelation Power', 'FontSize', 18);
end

% estimate dccFreq and check to see if there's a problem
interval = t(dccInd) - t(1);
dccFreq = 1.0 / interval;
if (dccFreq < 0.5 || dccFreq > 2.0) && ~(dccCorr > 0.3)
  % signal to other routines that there is a problem:
  dccCorr = NaN;
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dccFreq, dccPower, f, power] = refineDccFreq(t, dV, dccInd, ...
                                                       dccCorr, debugClean)
if dccInd <= 2  %Don't trust these results
  questionable = true;
  dccFreq = 1.0;  %Good a guess as any...
  dccRange = 2.0;  %Set a wide range
else
  Interval = t(dccInd) - t(1);
  dccFreq = 1.0 / Interval;
  if (dccFreq < 0.5 || dccFreq > 2.0) && ~(dccCorr > 0.3)
    dccFreq = 1.0;
    dccRange = 2.0;
    questionable = true;
  else
    questionable = false;
    dccRange = 1.25;
  end
end

if debugClean
  fprintf('Prelim:  Interval = %g, dccFreq = %g\n', Interval, dccFreq)
end

fSample = 1.0 / (t(2) - t(1));
try
  [power, f] = pmtm(dV, 7/2, [], fSample);
catch %#ok<CTCH>
  dV = dV(1:500000);
  [power, f] = pmtm(dV, 7/2, [], fSample);
end

ind = find(f > dccFreq / dccRange & f < dccFreq * dccRange);

[dccPower, ind2] = max(power(ind));
dccFreq = f(ind(ind2));

if isempty(dccFreq) || questionable
  disp('Warning:  QUESTIONABLE results in CleanDcc.m')
  dccFreq = Inf;
end

if debugClean
  fprintf('DCC Freq = %g kHz\n', dccFreq)
  h = NamedFigure('CleanDCC Power');
  set(h, 'WindowStyle', 'docked');
  plot(f, power);
  hold on
  plot(dccFreq, dccPower, 'ro')
  hold off
  xlabel('Frequency (kHz)', 'FontSize', 18)
  ylabel('Power', 'FontSize', 18)
  title('Clean DCC Power', 'FontSize', 18)
end
return