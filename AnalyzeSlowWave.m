function SlowWave = AnalyzeSlowWave(t, v, spikes, varargin)
% set the default options
defaultOptions = { ...
  'plotSubject', false, ...
  'doScalogram', false, ...
  'doSpectrum', true, ...
  'debugRefine', false, ...
  'debugPlots', false ...
};
% get the options overrides from varargin
options = GetOptions(defaultOptions, varargin);

if options.debugPlots && needPlot(options)
  options.plotGetSlowWave = options.plotSubject;
else
  options.plotGetSlowWave = false;
end

t = t / 1000;  %convert to seconds
DeltaT = t(2) - t(1);
SampleLen = t(end) - t(1);

[VWave, BackInds, SpikeAmp] = RemoveSpikes(v, spikes);

SlowWave = getPhases(VWave, DeltaT, BackInds, SpikeAmp, options);

SlowWave = CheckRealSlowWave(SlowWave, spikes, VWave, t);

if options.doSpectrum
  if options.doScalogram
    [spectrumInfo.Freq, spectrumInfo.Amplitude, spectrumInfo.Phase] ...
	    = Scalogram(v, t, options.plotSubject);
    spectrumInfo.AvgPower = GetSpectrumPower(spectrumInfo, SlowWave, [], ...
					 plotSubject, 'Average Power');
    SlowInds = GetSlowInds(spikes, t);
    if ~isempty(SlowInds)
      spectrumInfo.SlowPower = GetSpectrumPower(spectrumInfo, SlowWave, SlowInds, ...
					    plotSubject, ...
					    'Slow-wave Power');
    end
  else
    nw = 0.5 * round(2.0 * SampleLen / 30.0);
    if nw < 1
      nw = 1;
    end
    NumTapers = floor(2 * nw - 1);
    MaxFreq = 30.0;
    if size(v, 2) > 1
      v = v';
    end
    [spectrumInfo.Freq, spectrumInfo.Power, ~, ~, ~, ~, spectrumInfo.PowerConf] = ...
       coh_mt(DeltaT, v - mean(v), nw, NumTapers, MaxFreq, 0.68, 2);
    
    if options.debugPlots && needPlot(options)
      PSlow = interp1(spectrumInfo.Freq, spectrumInfo.Power, SlowWave.Freq);
      titleStr = makeTitle('Spectrum', options);
      h = NamedFigure(titleStr);
      set(h, 'WindowStyle', 'docked');
      clf; hold on
      %plot(Spectrum.Freq, sqrt(Spectrum.PowerConf(:,:,1)), 'r-');
      %plot(Spectrum.Freq, sqrt(Spectrum.PowerConf(:,:,2)), 'b-');
      plot(spectrumInfo.Freq, sqrt(spectrumInfo.Power), 'k-')
      plot(SlowWave.Freq, sqrt(PSlow), 'go', 'MarkerFaceColor', 'g');
      hold off
      title(RealUnderscores(titleStr), 'FontSize', 18);
      xlabel('Frequency (Hz)', 'FontSize', 18);
      ylabel('Power', 'FontSize', 18);
    end
  end
  SlowWave.Spectrum = spectrumInfo;
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vWave, n2List, spikeAmp] = RemoveSpikes(V, spikes)
vWave = V;
spikeAmp = 0;
n2List = [];
if isempty(spikes)  %Spikes weren't analyzed
  return
end
numSpikes = length(spikes.n1List);
if numSpikes == 0  %Spikes weren't found
  return
end

%n1List = Spikes.n1List;
%n2List = Spikes.n2List;
n1List = spikes.preMaxCurve.ind;
%n2List = spikes.postMaxCurve.ind;
%n1List = spikes.preMinV.ind;
n2List = spikes.postMinV.ind;

for SNum = 1:numSpikes
  n1 = n1List(SNum);
  n2 = n2List(SNum);
  vWave(n1:n2) = interp1([n1, n2], [V(n1), V(n2)], n1:n2);
  spikeAmp = spikeAmp + abs(V(n1) - V(n2));
end
spikeAmp = spikeAmp / numSpikes;
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SlowWave = getPhases(VWave, DeltaT, BackInds, SpikeAmp, options)

SigFact = 0.5;  %avg. amplitude must be larger than SigFact * SpikeAmp
SlowWave.Phases = [];
SlowWave.Amplitudes = [];
SlowWave.MinInds = [];

%First low-pass filter VWave (below 50 Hz) to eliminate noise
WFilt = (2 * DeltaT) * 50.0;
[B,A] = butter(2, WFilt, 'low');
%[B,A] = besself(1, WFilt, 'low');
VWave = filtfilt(B, A, VWave);

[SlowWave.Freq, SlowWave.Sigma, SlowWave.Corr] = ...
    GetSlowWave(DeltaT, VWave, options.plotGetSlowWave);
Freq = SlowWave.Freq;
if Freq <= 0
  plotWaveAnalysis(VWave, [], [], DeltaT, options);
  return
end

IndPeriod = round(1.0 / (DeltaT * Freq));
Phases = zeros(size(VWave));

NumWave = length(VWave);
LookFact = 0.1 + 0.4 * SlowWave.Corr;
LookAhead = round(LookFact * IndPeriod);

%Get a list of points smaller/larger than any point within distance LookAhead
[VWaveSort, SortInd] = sort(VWave);
MinInds = SortInd(1);
MaxInds = SortInd(end);
NumV = length(VWave);
m = NumV;
for n=2:NumV
  m = m - 1;
  Diff = min(abs(MinInds - SortInd(n)));
  if Diff > LookAhead
    MinInds = [MinInds, SortInd(n)];
  end
  Diff = min(abs(MaxInds - SortInd(m)));
  if Diff > LookAhead
    MaxInds = [MaxInds, SortInd(m)];
  end
end

%Remove minima that are too close to the back ends of spikes
if ~isempty(BackInds)
  n = 1;
  while n <= length(MinInds)
    m = MinInds(n);
    n2 = find(BackInds <= m, 1, 'last');
    if(isempty(n2))
      n = n + 1;
      continue
    end
    n2 = BackInds(n2);
    if m - n2 <= 1 %Remove this MinInd
      MinInds = [MinInds(1:(n-1)), MinInds((n+1):end)];
    else
      n = n + 1;
    end
  end
end

%Refine the list to only include local minima/maxima
n = 2;
while n <= length(MinInds)
  Test = MinInds(n);
  I1 = max([1, Test - LookAhead]);
  I2 = min([NumV, Test + LookAhead]);
  LocalMin = min(VWave(I1:I2));
  if(LocalMin < VWave(Test))
    MinInds = [MinInds(1:(n-1)), MinInds((n+1):end)];
  else
    n = n + 1;
  end
end
MinInds = sort(MinInds);
n = 2;
while n <= length(MaxInds)
  Test = MaxInds(n);
  I1 = max([1, Test - LookAhead]);
  I2 = min([NumV, Test + LookAhead]);
  LocalMax = max(VWave(I1:I2));
  if(LocalMax > VWave(Test))
    MaxInds = [MaxInds(1:(n-1)), MaxInds((n+1):end)];
  else
    n = n + 1;
  end
end
MaxInds = sort(MaxInds);
LocalMins = VWave(MinInds);
LocalMaxes = VWave(MaxInds);

if length(MinInds) < 3 || length(MaxInds) < 3
  plotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, options);
  return
end

%now make sure they alternate
if(MinInds(1) < MaxInds(1))
  Current = -1;
else
  Current = 1;
end
nMin = 1;
nMax = 1;
done = false;
while ~done
  if Current > 0
    if nMax == length(MaxInds)
      done = true;
      %prevent several mins in a row at the end:
      if nMin < length(MinInds)
        [TempMin, TempInd] = min(VWave((MaxInds(end)+1):end));
        TempInd = TempInd + MaxInds(end);
        nMin = nMin - 1;
        MinInds = [MinInds(1:nMin), TempInd];
        LocalMins = [LocalMins(1:nMin), TempMin];
      end
      break
    end
    if MaxInds(nMax + 1) < MinInds(nMin)
      %Two maxes in a row.  Remove smaller
      if LocalMaxes(nMax) > LocalMaxes(nMax + 1)
        nDel = nMax + 1;
      else
        nDel = nMax;
      end
      MaxInds = [MaxInds(1:(nDel - 1)), MaxInds((nDel + 1):end)];
      LocalMaxes = [LocalMaxes(1:(nDel - 1)), LocalMaxes((nDel + 1):end)];
    else
      nMax = nMax + 1;
      Current = -1;
    end
  else
    if nMin == length(MinInds)
      done = true;
      %prevent several maxes in a row at the end:
      if nMax < length(MaxInds)
        [TempMax, TempInd] = max(VWave((MinInds(end)+1):end));
        TempInd = TempInd + MinInds(end);
        nMax = nMax - 1;
        MaxInds = [MaxInds(1:nMax), TempInd];
        LocalMaxes = [LocalMaxes(1:nMax), TempMax];
      end
      break
    end
    if MinInds(nMin + 1) < MaxInds(nMax)
      %Two mins in a row.  Remove bigger
      if LocalMins(nMin) < LocalMins(nMin + 1)
        nDel = nMin + 1;
      else
        nDel = nMin;
      end
      MinInds = [MinInds(1:(nDel - 1)), MinInds((nDel + 1):end)];
      LocalMins = [LocalMins(1:(nDel - 1)), LocalMins((nDel + 1):end)];
    else
      nMin = nMin + 1;
      Current = 1;
    end
  end
end

if length(MinInds) < 3 || length(MaxInds) < 3
  plotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, options);
  return
end

%ensure that each min is a local min between maxes, and vice versa
changed = true;
while changed
  changed = false;
  if MinInds(1) < MaxInds(1)
    Current = 1;
  else
    Current = -1;
  end
  nMin = 1;
  nMax = 1;
  done = false;
  while ~done
    if Current > 0
      if nMin == length(MinInds)
        done = true;
        I2 = NumV;
      else
        I2 = MinInds(nMin + 1) - 1;
      end
      I1 = MinInds(nMin) + 1;
      [Val, Ind] = max(VWave(I1:I2));
      if Val > LocalMaxes(nMax)
        %max out of place
        changed = true;
        MaxInds(nMax) = Ind + I1 - 1;
        LocalMaxes(nMax) = Val;
      end
      nMin = nMin + 1;
      Current = -1;
    else
      if nMax == length(MaxInds)
        done = true;
        I2 = NumV;
      else
        I2 = MaxInds(nMax + 1) - 1;
      end
      I1 = MaxInds(nMax) + 1;
      [Val, Ind] = min(VWave(I1:I2));
      if Val < LocalMins(nMin)
        %max out of place
        changed = true;
        MinInds(nMin) = Ind + I1 - 1;
        LocalMins(nMin) = Val;
      end
      nMax = nMax + 1;
      Current = 1;
    end
  end
end

if length(MinInds) < 3 || length(MaxInds) < 3
  plotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, options);
  return
end

%Remove min/max pairs until the slow wave frequency seems to match
DeltaFreq = 1.0 / (DeltaT * (MinInds(end) - MinInds(1)));
Freq = (length(MinInds) - 1) * DeltaFreq;
DiffMaxes = GetDiffMeasure(LocalMaxes, 1, options.debugRefine);
DiffMins = GetDiffMeasure(LocalMins, -1, options.debugRefine);
if options.debugRefine
  fprintf('CurrentFreq: %g, Actual: %g\n', Freq, SlowWave.Freq);
end

while Freq - SlowWave.Freq > DeltaFreq || DiffMaxes >= 3 || DiffMins >= 3
  if DiffMaxes > DiffMins
    [~, RMaxInd] = min(LocalMaxes);
    IndMax = MaxInds(RMaxInd);
    Ind2 = find(MinInds > IndMax, 1);
    if isempty(Ind2)
      RMinInd = length(LocalMins);
    else
      Ind1 = Ind2 - 1;
      if Ind1 == 0 || LocalMins(Ind1) < LocalMins(Ind2)
        RMinInd = Ind2;
      else
        RMinInd = Ind1;
      end
    end
  else
    [~, RMinInd] = max(LocalMins);
    IndMin = MinInds(RMinInd);
    Ind2 = find(MaxInds > IndMin, 1);
    if isempty(Ind2)
      RMaxInd = length(LocalMaxes);
    else
      Ind1 = Ind2 - 1;
      if Ind1 == 0 || LocalMaxes(Ind1) > LocalMaxes(Ind2)
        RMaxInd = Ind2;
      else
        RMaxInd = Ind1;
      end    
    end
  end
  
  MinInds = [MinInds(1:(RMinInd-1)), MinInds((RMinInd+1):end)];
  LocalMins = [LocalMins(1:(RMinInd-1)), LocalMins((RMinInd+1):end)];
  MaxInds = [MaxInds(1:(RMaxInd-1)), MaxInds((RMaxInd+1):end)];
  LocalMaxes = [LocalMaxes(1:(RMaxInd-1)), LocalMaxes((RMaxInd+1):end)];
  
  if length(MinInds) < 3 || length(MaxInds) < 3
    plotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, options);
    return
  end

  DeltaFreq = 1.0 / (DeltaT * (MinInds(end) - MinInds(1)));
  Freq = (length(MinInds) - 1) * DeltaFreq;
  if options.debugRefine
    fprintf('CurrentFreq: %g, Actual: %g\n', Freq, SlowWave.Freq);
  end
  DiffMaxes = GetDiffMeasure(LocalMaxes, 1, options.debugRefine);
  DiffMins = GetDiffMeasure(LocalMins, -1, options.debugRefine);
end

%MinInds and MaxInds are calculated
%  now calculate Phases:
if MinInds(1) < MaxInds(1)
  Current = -1;
  CurrentMax = LocalMaxes(1);
  nMin = 1;
  nMax = 0;
else
  Current = 1;
  CurrentMin = LocalMins(1);
  nMin = 0;
  nMax = 1;
end
done = false;
IndStart = 1;
Amplitudes = [];
while ~done
  if Current > 0  %increasing
    if nMax > length(MaxInds)
      IndStop = NumV;
      done = true;
    else
      CurrentMax = LocalMaxes(nMax);
      IndStop = MaxInds(nMax);
    end
    Ind = IndStart:IndStop;
    IndStart = IndStop + 1;
    Slope = 2.0 / (CurrentMax - CurrentMin);
    Offset = 1.0 + CurrentMin * Slope;
    Phases(Ind) = 0.5 * pi + real(asin(Slope * VWave(Ind) - Offset));
    nMin = nMin + 1;
    Current = -1;
  else   %decreasing
    if nMin > length(MinInds)
      IndStop = NumV;
      done = true;
    else
      CurrentMin = LocalMins(nMin);
      IndStop = MinInds(nMin);
    end
    Ind = IndStart:IndStop;
    IndStart = IndStop + 1;
    Slope = 2.0 / (CurrentMax - CurrentMin);
    Offset = 1.0 + CurrentMin * Slope;
    Phases(Ind) = 1.5 * pi - real(asin(Slope * VWave(Ind) - Offset));
    nMax = nMax + 1;
    Current = 1;
    if nMax <= length(MaxInds)
      Amplitudes = [Amplitudes, ...
                    0.5 * (LocalMaxes(nMax) + CurrentMax) - CurrentMin];
    end
  end
end

%fprintf('MeanAmp: %g, SpikeAmp: %g, CutOff:  %g\n', ...
%	mean(Amplitudes), SpikeAmp, SigFact * SpikeAmp)
if mean(Amplitudes) < SigFact * SpikeAmp
  plotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, options);
  return
end
SlowWave.Phases = Phases;
SlowWave.Amplitudes = Amplitudes;
SlowWave.MinInds = MinInds;
plotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, options);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotWaveAnalysis(VWave, MinInds, MaxInds, DeltaT, options)
IntensePlot = false;
if options.debugPlots && needPlot(options)
  titleStr = 'Slow-wave Phase';
  if ischar(options.plotSubject)
    titleStr = [options.plotSubject, ' - ', titleStr];
  end
  h = NamedFigure(titleStr);
  set(h, 'WindowStyle', 'docked');
  clf
  hold on
  if IntensePlot
    t = 0;
    for n = 1:length(VWave)
      Blue = 0.5*(1.0 + cos(Phases(n)));
      Red = 1.0 - Blue;
      plot(t, VWave(n), '.', 'MarkerEdgeColor', [Red, 0, Blue], ...
	   'MarkerSize', 1)
      t = t + DeltaT;
    end
  else
    plot(DeltaT * (0:(length(VWave)-1)), VWave, 'b-');
  end
  if length(MaxInds) >= 3 && length(MinInds) >= 3
    plot(MaxInds * DeltaT, VWave(MaxInds), 'go');
    plot(MinInds * DeltaT, VWave(MinInds), 'gx');
  end
  hold off
  title(RealUnderscores(titleStr), 'FontSize', 18);
  xlabel('Time (s)', 'FontSize', 18);
  ylabel('Voltage (mV)', 'FontSize', 18);
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DiffMeasure = GetDiffMeasure(X, Type, debugRefine)
%  Type is +1 for Maxes, -1 for Mins
d = pdist(X');
Z = linkage(d);
c = cluster(Z, 'maxclust', 2);

Ind1 = find(c == 1);
Mean1 = mean(X(Ind1));
Std1 = std(X(Ind1));
Ind2 = find(c == 2);
Mean2 = mean(X(Ind2));
Std2 = std(X(Ind2));

%look for abnormally small Maxes, abnormally large Mins
if Mean1 * Type > Mean2 * Type
  L1 = length(Ind1);
  L2 = length(Ind2);
else
  L1 = length(Ind2);
  L2 = length(Ind1);
end
if L1 > L2 || L1 >= 3  %could be abnormal in way we are concerned
  DiffScale = 2 * max(Std1, Std2);
  if(DiffScale < 1)
    DiffScale = 1;
  end
  DiffMeasure = abs(Mean2 - Mean1) / DiffScale;
else   %abnormally large Maxes or abnormally small Mins, don't care:
  DiffMeasure = 0;
end

if debugRefine
  fprintf('M1: %g, S1: %g, M2: %g, S2: %g, Diff: %g\n', ...
          Mean1, Std1, Mean2, Std2, DiffMeasure);
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MinInds, MaxInds] = CheckNew(VWave, MinInds, MaxInds, Type, ...
				       IndMin, IndMax)
if nargin < 6
  IndMin = length(MinInds);
  IndMax = length(MaxInds);
end
if strcmp(Type, 'min')
  if(IndMin < 2)
    return
  end
  Ind = MinInds(IndMin - 1):MinInds(IndMin);
  [Val, NewMaxInd] = max(VWave(Ind));
  NewMaxInd = NewMaxInd + Ind(1) - 1;
  if(MaxInds(IndMax) ~= NewMaxInd)
    MaxInds(IndMax) = NewMaxInd;
    [MinInds, MaxInds] = CheckNew(VWave, MinInds, MaxInds, 'max', ...
				  IndMin - 1, IndMax);
  end
else
  if IndMax < 2
    return
  end
  Ind = MaxInds(IndMax - 1):MaxInds(IndMax);
  [Val, NewMinInd] = min(VWave(Ind));
  NewMinInd = NewMinInd + Ind(1) - 1;
  if(MinInds(IndMin) ~= NewMinInd)
    MinInds(IndMin) = NewMinInd;
    [MinInds, MaxInds] = CheckNew(VWave, MinInds, MaxInds, 'min', ...
				  IndMin, IndMax - 1);
  end  
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SlowWave = CheckRealSlowWave(SlowWave, Spikes, VWave, t)
%First check to see if the slow-wave is significantly larger than
%  background noise

if(SlowWave.Freq <= 0)
  SlowWave.NumBurst = 0;
  SlowWave.NumOnlySpike = 0;
  SlowWave.NumNoSpike = 0;
  SlowWave.SpikesPerWave = [];
  return
end

%high-pass filter the waveform
SmoothTime = .1 * (1.0 / SlowWave.Freq);
n = ceil(SmoothTime / (t(2) - t(1)));
if(n < 8)
  n = 8;
end
[B,A] = butter(2, 2 / n, 'high');  %2/n because Nyquist rate is 1/2
%[B,A] = besself(4, 2 / n, 'high');  %2/n because Nyquist rate is 1/2
VFilt = filtfilt(B, A, VWave);

VFilt = sort(abs(VFilt));
OneSigma = 0.682689492;
TwoSigma = 0.954499736;
Fuzz = 2 * VFilt(round(OneSigma * length(VFilt)));  %Factor of two
                                                    %makes two-sided

MedianAmp = median(SlowWave.Amplitudes);
SlowWave.Sigma = MedianAmp / Fuzz;
if SlowWave.Sigma < 2
  SlowWave.Amplitudes = [];
  SlowWave.MinInds = [];
  SlowWave.NumBurst = 0;
  SlowWave.NumOnlySpike = 0;
  SlowWave.NumNoSpike = 0;
  SlowWave.SpikesPerWave = [];
  return
end
  
%Now check to see if spikes are the cause of the slow-wave.
% If they are, the phase should be advanced (near trough of slow-wave)
% at the end of the spike.
if isempty(Spikes)
  NumSpikes = 0;
else
  n1List = Spikes.n1List;
  NumSpikes = length(n1List);
end
NumWaves = length(SlowWave.MinInds) - 1;
if(NumSpikes == 0)
  SlowWave.NumBurst = 0;
  SlowWave.NumOnlySpike = 0;
  SlowWave.NumNoSpike = NumWaves;
  SlowWave.SpikesPerWave = zeros(NumWaves, 1);
  return
end
n2List = Spikes.n2List;
spikeTimes = Spikes.times;
Phases = SlowWave.Phases;
MinInds = SlowWave.MinInds;

PhaseCutoff = 1.5 * pi;

NumBurst = 0;  %Increment if the phase is not advanced
NumOnlySpike = 0;  %Increment if the phase is advanced
NumNoSpike = 0;  %Increment if there is no spike at all
SpikesPerWave = zeros(NumWaves, 1);
for WaveNum = 1:NumWaves
  StartInd = MinInds(WaveNum);
  StopInd = MinInds(WaveNum + 1);
  startT = 1000.0 * t(StartInd); stopT = 1000.0 * t(StopInd);
  SpikeInds = find(startT < spikeTimes & stopT > spikeTimes);
  Num = length(SpikeInds);
  SpikesPerWave(WaveNum) = Num;
  if(Num == 0)
    NumNoSpike = NumNoSpike + 1;
  elseif(Num == 1)
    n2 = n2List(SpikeInds(end));
    Theta_n2 = Phases(n2);
    MinPhase = min(Phases(n2:StopInd-1));
    %fprintf('p_n2: %g, MinP: %g\n', Theta_n2, MinPhase))
    if(MinPhase > PhaseCutoff)
      NumOnlySpike = NumOnlySpike + 1;
    else
      NumBurst = NumBurst + 1;
    end
  else
    NumBurst = NumBurst + 1;
  end
end

SlowWave.NumBurst = NumBurst;
SlowWave.NumOnlySpike = NumOnlySpike;
SlowWave.NumNoSpike = NumNoSpike;
SlowWave.SpikesPerWave = SpikesPerWave;
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Power = GetSpectrumPower(Spectrum, SlowWave, Indices, ...
				  PlotSubject, TitleStr)
if isempty(Indices)
  Power = mean(Spectrum.Amplitude.^2, 2);
else
  Power = mean(Spectrum.Amplitude(:,Indices).^2, 2);
end

if needPlot(PlotSubject)
  PSlow = interp1(Spectrum.Freq, Power, SlowWave.Freq);
  if ischar(PlotSubject) && ~isempty(PlotSubject)
    TitleStr = [PlotSubject, ': ', TitleStr];
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  hold off
  loglog(Spectrum.Freq, Power)
  hold on
  plot(SlowWave.Freq, PSlow, 'ro');
  hold off
  title(RealUnderscores(TitleStr), 'FontSize', 18);
  xlabel('Frequency (Hz)', 'FontSize', 18);
  ylabel('Power', 'FontSize', 18);
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotVar = needPlot(options)
plotVar = ischar(options.plotSubject) || options.plotSubject;
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function titleStr = makeTitle(baseTitle, options)
titleStr = baseTitle;
if ischar(options.plotSubject)
  titleStr = [options.plotSubject, ' - ', titleStr];
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SlowInds = GetSlowInds(Spikes, t)
n1List = Spikes.n1List;
n2List = Spikes.n2List;
NumT = length(t);
if(isempty(n1List))
  SlowInds = 1:NumT;
  return
end

if(n1List(1) > 1)
  SlowInds = 1:(n1List(1)-1);
else
  SlowInds = [];
end
SlowInds = [];

for n=2:length(n1List)
  StartInd = n2List(n-1) + 1;
  StopInd = n1List(n) - 1;
  Mid = round(0.5* (StartInd + StopInd));
  HalfRange = round(0.2 * 0.5 * (StopInd - StartInd));
  StartInd = Mid - HalfRange;
  StopInd = Mid + HalfRange;
  
  if(StopInd >= StartInd)
    SlowInds = [SlowInds, StartInd:StopInd];
  end
end

%if(n2List(end) < NumT)
%  SlowInds = [SlowInds, (n2List(end)+1):NumT];
%end

return