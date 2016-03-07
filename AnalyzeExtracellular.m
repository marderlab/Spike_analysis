function structOut = AnalyzeExtracellular(t, v, varargin)
% structOut = AnalyzeExtracellular(t, v, plotSubject)
% Analyzes a single extracellular waveform, returning
% information about bursting

callstack = dbstack;
if length(callstack) == 1  % not called by another function
  tic
end

if length(t) ~= length(v)
  error('Time and Voltage arrays have different length!')
elseif size(t, 1) > 1
  t = t';
  v = v';
end

% set the default options
defaultOptions = { ...
  'plotSubject', false, ...
  'lowCutoff', NaN, ...
  'highCutoff', NaN, ...
  'bracketWidth', 15.0, ...
  'pFalseSpike', 0.05, ...
  'debugPlots', false ...
};
% get the options overrides from varargin
options = GetOptions(defaultOptions, varargin);

spike = getSpike(t, v, options);
burst = getBurst(t, v, spike, options);

structOut.spike = spike;
structOut.burst = burst;

if needPlot(options)
  plotAnalyze(t, v, structOut, options)
end

if length(callstack) == 1  % not called by another function
  Toc
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spike = getSpike(t, v, options)
maxTimeWidth = 15.0; % ms

numV = length(v);
dT = (t(end) - t(1)) / (numV - 1);
indWidth = ceil(maxTimeWidth / dT);
%[lowCutoff, highCutoff] = getCutoffs(v);
[lowCutoff, highCutoff] = getAutoCutoffs(dT, v, options);

n1List = [];
n2List = [];

n1 = 1;
while n1 < numV
  % try to bracket a spike between n1 and n2, as a large positive
  % deflection at n1, followed by a large negative deflection at n2
  if v(n1) < highCutoff
    n1 = n1 + 1;
    continue
  end
  
  n2 = n1 + 1;
  n2Stop = min(n1 + indWidth, numV);
  while n2 <= n2Stop
    if v(n2) < lowCutoff
      break
    elseif v(n2) > highCutoff
      n1 = n2;
      n2Stop = min(n1 + indWidth, numV);
    end
    
    n2 = n2 + 1;
  end
  
  if n2 > n2Stop
    n1 = n2;
    continue
  else
    n1 = n2 + 1;
  end
  
  % spike exists between n1 and n2
  n1List = [n1List, n1];
  n2List = [n2List, n2];
end

spike.n1List = n1List;
spike.n2List = n2List;
spike.times = 0.5 * (t(n1List) + t(n2List));
numSpikes = length(n1List);
if numSpikes == 0
  spike.frequency = 0;
  spike.ISI = [];
elseif numSpikes == 1
  spike.frequency = 1.0 / (t(end) - t(1) + dT);
  spike.ISI = [];
else
  spike.frequency = (numSpikes - 1) / ...
      (spike.times(end) - spike.times(1));
  spike.ISI = diff(spike.times);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lowCutoff, highCutoff] = getAutoCutoffs(dT, deriv, options)
% Get cutoffs for significant spiking

% first compute how rare/extreme derivatives cutoffs need to be to set a
% false detection probability

% number of trace points in a bracketed spike
nBracket = options.bracketWidth / dT; % number of points in a brac
len = length(deriv);
logOdds = 4 * log(1 - options.pFalseSpike) / len / nBracket;

% this is how rare a derivative has to be (either positive or negative) to
% achieve the given false-detection probability
minRareness = sqrt(-logOdds);

% compute approximate 1-sigma levels for positive and negative derivatives,
% based on presumably nearly-gaussian small derivatives
sortDeriv = sort(deriv);
medianInd = round(0.5 * length(sortDeriv));
medianDV = sortDeriv(medianInd);
sigmaFact = 0.5 * erf(1.0/sqrt(2));
sigmaPosInd = round((0.5 + sigmaFact) * length(sortDeriv));
sigmaPos = sortDeriv(sigmaPosInd) - medianDV;
sigmaNegInd = round((0.5 - sigmaFact) * length(sortDeriv));
sigmaNeg = medianDV - sortDeriv(sigmaNegInd);

% estimate number of sigma needed to achieve minRareness
numSigma = sqrt(2) * erfcinv(2 * minRareness);

% compute cutoffs by moving them appropriate number of sigma away from
% median
lowCutoff = medianDV - numSigma * sigmaNeg;
highCutoff = medianDV + numSigma * sigmaPos;

if options.debugPlots
  titleStr = makeTitle('Spike Thresholds', options);
  
  h = NamedFigure(titleStr);
  set(h, 'WindowStyle', 'docked')
  clf
  [n, x] = hist(sortDeriv, 1000);
  n = n ./ max(n);
  bar(x, n);
  hold on
  plot([lowCutoff, lowCutoff], [0, 1], 'r')
  plot([highCutoff, highCutoff], [0, 1], 'g')
  hold off
  xlabel('Derivative (mV/ms)')
  ylabel('Relative Frequency')
  title(RealUnderscores(titleStr))
  legend('Derivatives', 'Low threshold', 'High threshold')
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lowCutoff, highCutoff] = getCutoffs(v)
% Get cutoffs for significant spiking
cutoffFrac = 0.999; % Rareness is setting the scale for "large" events
discount = 0.50;  % Since rareness is significant, we want to allow
                 % things to be somewhat smaller than the rarest
                 % events and still be significant
meanV = mean(v);
		 
vNeg = sort(v(find(v < meanV)), 'descend');
vPos = sort(v(find(v > meanV)));

indNeg = round(cutoffFrac * length(vNeg));
indPos = round(cutoffFrac * length(vPos));

lowCutoff = meanV + discount * (vNeg(indNeg) - meanV);
highCutoff = meanV + discount * (vPos(indPos) - meanV);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function burst = getBurst(t, v, spike, options)
maxSpacing = 200; %ms

ind = find(spike.ISI > maxSpacing);

n1List = spike.n1List(1 + ind(1:(end-1)));
n2List = spike.n2List(ind(2:end));

% check if first and last spike are far enough from edges
%  to qualify as bursts
if t(spike.n1List(1)) - t(1) > maxSpacing
  n1List = [spike.n1List(1), n1List];
  n2List = [spike.n2List(ind(1)), n2List];
end
if t(end) - spike.n2List(end) > maxSpacing
  n1List = [n1List, spike.n1List(1 + ind(end))];
  n2List = [n2List, spike.n2List(end)];
end

burst.n1List = n1List;
burst.n2List = n2List;
burst.startTimes = t(n1List);
burst.stopTimes = t(n2List);

burst.durations = structifyList(burst.stopTimes - burst.startTimes);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outStruct = structifyList(inList)
outStruct.List = inList;
goodInd = find(isfinite(inList));
if length(goodInd) > 1
  inList = inList(goodInd);
  outStruct.Mean = mean(inList);
  outStruct.StdDev = std(inList);
  outStruct.Variance = outStruct.StdDev^2;
  outStruct.CoefOfVar = outStruct.StdDev / outStruct.Mean;
else
  if length(goodInd) == 1
    outStruct.Mean = inList(goodInd);
  else
    outStruct.Mean = 0;
  end
  outStruct.StdDev = 0;
  outStruct.Variance = 0;
  outStruct.CoefOfVar = 0;
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotVar = needPlot(options)
plotVar = ischar(options.plotSubject) || options.plotSubject;
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAnalyze(t, v, structOut, options)
if ischar(options.plotSubject)
  titleStr = options.plotSubject;
else
  titleStr = 'Analyze Extracellular';
end

maxV = max(v);
minV = min(v);
bottom = minV;
top = maxV;

h = NamedFigure(titleStr);
set(h, 'WindowStyle', 'docked')
clf
whitebg(h, 'k');
hold on
% plot bursts
numBursts = length(structOut.burst.startTimes);
for n = 1:numBursts
  tLow = 0.001 * structOut.burst.startTimes(n);
  tHigh = 0.001 * structOut.burst.stopTimes(n);
  fill([tLow, tHigh, tHigh, tLow], [bottom, bottom, top, top], 'b');
end

% plot spikes
numSpikes = length(structOut.spike.times);
for n = 1:numSpikes
  t_n = 0.001 * structOut.spike.times(n);
  plot([t_n, t_n], [bottom, top], 'r-')
end

% plot the waveform
plot(0.001 * t, v, 'w-');

hold off

title(titleStr, 'FontSize', 18)
xlabel('Time (s)', 'FontSize', 18)
ylabel('Voltage (arbitrary)', 'FontSize', 18)
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