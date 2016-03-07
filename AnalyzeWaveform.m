function structOut = AnalyzeWaveform(t, v, varargin)
% structOut = AnalyzeWaveform(t, v, plotSubject)
% Analyzes a single voltage waveform, looking for spikes
%    and bursts, and calculating relevant frequencies.
%
%  INPUT PARAMETERS:
%   -t is time in ms
%   -v is voltage in mV
%    OPTIONAL:
%     -plotSubject should be set to true[false] to produce[suppress]
%       plots of waveforms/analysis.  Alternatively, it can be set
%       to a string to aid it titling plots (e.g. 'Exp #71')
%       plotSubject defaults to false
%  OUTPUT PARAMETERS:
%   -structOut.spikes:  structure with spike information
%      -spikes.Freq is overall spiking frequency (in Hz)
%      -spikes.Times is a plain list of spike times (in ms)
%      -spikes.Intervals is a list of interspike intervals (in ms)
%      -spikes.Frequencies is a list of instantaneous frequencies (in Hz)
%      Shape information structures (should be self-descriptive)
%      -spikes.MaxV, Spike.MaxDeriv, Spike.MinDeriv, Spike.PreMinV,
%       spikes.PostMinV, Spike.PreMaxCurve, Spike.PostMaxCurve
%           Each contains a list of times/voltage points, and if relevant
%           another quantity (such as K for curvatures)
%   -structOut.slowWave:  structure with slow-wave information
%      -slowWave.Freq: frequency of the dominant slow-wave
%       component (in Hz)
%      -slowWave.Sigma: (very crude) measure of the importance
%       of the slow-wave frequency in the power spectrum
%      -slowWave.Corr: autocorrelation at slow-wave period
%      -slowWave.Spectrum: structure with spectrum information
%        -Spectrum.Freq:  list of analyzed frequencies
%        -Spectrum.Power: length NumFreq list of average powers of
%           waveform with spikes removed.
%   -structOut.bursts:  structure with burst information
%      -bursts.Freq is burst frequency (in Hz)
%      -bursts.SpikeFreq is within-burst spike frequency (in Hz)
%      -bursts.DutyCycle is the average burst duration/period
%      -bursts.Times is a plain list of burst times (in ms)
%      -bursts.Durations is a list of burst durations (in ms)
%      -bursts.numSpikes is a list of spikes per burst
%      -bursts.SpikeFrequencies is a list of spike frequencies (in Hz)
%      -bursts.InterBurstIntervals is a list of inter-burst
%       intervals (in ms)
%   -structOut.medianV:  the median of the voltage trace.  If the
%      cell is silent, it should be the resting potential,
%      otherwise, who knows...
%
%List structures usually will have a Name.List element, as well as
%  Name.Mean, Name.StdDev, Name.Variance, Name.CoefOfVar
%  (a few are just plain lists)
%If a feature is not detected, relevant frequencies are set to
%  zero, and relevant lists are empty
%
%NOTE for future:  would benefit enormously by changing to .mex
callstack = dbstack;
if length(callstack) == 1  % not called by another function
  tic
end

if nargin < 2
  help AnalyzeWaveform
  error('Invalid number of arguments')
end
if length(t) ~= length(v)
  if length(t) == 1
    dt = t;
    t = 0:dt:(dt * (length(v) - 1));
  else
    error('Time and Voltage arrays have different length!')
  end
end
if size(t, 1) > 1
  t = t';
end
if size(v,2) ~= size(t,2)
  v = v';
end

% set the default options
defaultOptions = { ...
  'plotSubject', false, ...
  'timesOnly', false, ...
  'firstOnly', false, ...
  'lowCutoff', NaN, ...
  'highCutoff', NaN, ...
  'bracketWidth', 5.0, ...
  'minCutoffDiff', 0.1, ...
  'pFalseSpike', 1.0e-4, ...
  'recursive', false, ...
  'removeOutliers', true, ...
  'findMinis', false, ...
  'debugPlots', false ...
};
% get the options overrides from varargin
options = GetOptions(defaultOptions, varargin, true);

spikes = GetSpikes(t, v, varargin{:});

slowWave = AnalyzeSlowWave(t, v, spikes, ...
  'plotSubject', options.plotSubject, 'debugPlots', options.debugPlots);
bursts = AnalyzeBurst(spikes, slowWave, t);
slowWave.Phases = [];  %Reduce storage demand

%structify (add info about mean, variance, etc) various lists
spikes.intervals = structifyList(spikes.intervals);
spikes.frequencies = structifyList(spikes.frequencies);

spikes.maxV.v = structifyList(spikes.maxV.v);
spikes.maxDeriv.v = structifyList(spikes.maxDeriv.v);
spikes.maxDeriv.dV = structifyList(spikes.maxDeriv.dV);
spikes.minDeriv.v = structifyList(spikes.minDeriv.v);
spikes.minDeriv.dV = structifyList(spikes.minDeriv.dV);
spikes.preMinV.v = structifyList(spikes.preMinV.v);
spikes.postMinV.v = structifyList(spikes.postMinV.v);
spikes.preMaxCurve.v = structifyList(spikes.preMaxCurve.v);
spikes.preMaxCurve.K = structifyList(spikes.preMaxCurve.K);
spikes.postMaxCurve.v = structifyList(spikes.postMaxCurve.v);
spikes.postMaxCurve.K = structifyList(spikes.postMaxCurve.K);

bursts.Durations = structifyList(bursts.Durations);
bursts.SpikesPerBurst = structifyList(bursts.SpikesPerBurst);
bursts.InterBurstIntervals = structifyList(bursts.InterBurstIntervals);
bursts.SpikeFrequencies = structifyList(bursts.SpikeFrequencies);

structOut.spikes = spikes;
structOut.slowWave = slowWave;
structOut.bursts = bursts;
structOut.medianV = median(v);

if needPlot(options)
  dT = t(2) - t(1);
  hSpikes = PlotGetSpikes(dT, v, spikes, options, bursts);
  
  % link relevant time axis together
  if options.debugPlots
    aSpikes = get(hSpikes, 'CurrentAxes');
    derivsTitle = makeTitle('dV/dT vs. t', options);
    aDerivs = get(findobj('name', derivsTitle),'CurrentAxes');
    aHandles = [aSpikes, aDerivs];
    linkaxes(aHandles, 'x');
  end
end

if length(callstack) == 1  % not called by another function
  Toc
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outStruct = structifyList(inList)
outStruct.list = inList;
goodInd = find(isfinite(inList));
if length(goodInd) > 1
  inList = inList(goodInd);
  outStruct.mean = mean(inList);
  outStruct.stdDev = std(inList);
  outStruct.variance = outStruct.stdDev^2;
  outStruct.coefOfVar = outStruct.stdDev / outStruct.mean;
else
  if length(goodInd) == 1
    outStruct.mean = inList(goodInd);
  else
    outStruct.mean = 0;
  end
  outStruct.stdDev = 0;
  outStruct.variance = 0;
  outStruct.coefOfVar = 0;
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotVar = needPlot(options)
plotVar = ischar(options.plotSubject) || options.plotSubject;
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function titleStr = makeTitle(titleBase, options)
% set the full title for a figure based on base title and plotSubject
if ischar(options.plotSubject)
  titleStr = [options.plotSubject, ': ', titleBase];
else
  titleStr = titleBase;
end
return