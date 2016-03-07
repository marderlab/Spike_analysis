function [Freq, varargout] = GetSlowWave(DeltaT, V0, V1, PlotSubject)
%[Freq, NumSigma, SlowC] = GetSlowWave(DeltaT, V0, V1, PlotSubject)
% calculates the slow wave frequency and sigma
%   INPUT PARAMETERS:
%    -DeltaT is time step (s)
%    -V0 is first voltage trace (arbitrary units)
%    OPTIONAL:
%     -V1 is second voltage trace (arbitrary units)
%     -PlotSubject  set to 1/true to plot.  Defaults to false
%   OUTPUT PARAMETERS:
%     -Freq is dominant slow wave frequency in Hz
%      OPTIONAL:
%      -NumSigma is z-score of slow wave frequency correlation,
%       compared to average
%      -SlowC is autocorrelation at slow wave period

%  METHODS write-up
% Slow wave detection.  Here we mean low-frequency membrane oscillations,
% agnostic as to origin.  We identified the frequency of these oscillations
% as the lowest-frequency peak in the autocorrelogram of the voltage
% waveform.  If this peak had an autocorrelation less than 0.05, we instead
% used the frequency corresponding to the greatest autocorrelation (this is
% uncommon and corresponds to very irregular activity). 

CCutoff = 0.05;
%First calculate the autocorrelations:
NumAutoCorr = round(length(V0) / 3);  %need to see at least 3 waves!
C = xcorr(zscore(V0), NumAutoCorr, 'unbiased');
C = C((NumAutoCorr+1):end);
f = (1.0 / DeltaT) ./ (1:length(C));

if nargin >= 3
  if length(V1) == length(V0)
    C1 = xcorr(zscore(V1), NumAutoCorr, 'unbiased');
    C1 = C1((NumAutoCorr+1):end);
    C = 0.5 * (C + C1);
    if nargin == 3
      PlotSubject = false;
    end
  elseif nargin == 3
    PlotSubject = V1;
  else
    error('Invalid value for V1')
  end
else
  PlotSubject = false;
end

%Find the slow wave peak.
[MaxInd, MaxC, IndStart, IndStop] = FindSlowWavePeak(C);
Freq = f(MaxInd);

if isempty(Freq)
  Freq = 0;
  MaxC = 0;
  IndStart = 1;
  IndStop = 1;
elseif MaxC < CCutoff
  TempStart = IndStart;
  IndStart = find(C < CCutoff, 1);
  IndStart = find(C((IndStart + 1):end) >= CCutoff, 1) + IndStart;
  if isempty(IndStart)
    IndStart = TempStart;
  else
    IndStop = find(C((IndStart+1):end) < CCutoff, 1) + IndStart;
    if isempty(IndStop)
      IndStop = length(C);
    end
    [MaxC, MaxInd] = max(C(IndStart:IndStop));
    MaxInd = MaxInd + IndStart - 1;
    Freq = f(MaxInd);
  end
end
if DoPlot(PlotSubject)
  TitleStr = 'SlowWave Correlation';
  if ischar(PlotSubject) && ~isempty(PlotSubject)
    TitleStr = [PlotSubject, ' - ', TitleStr];
  else
    
  end
  h = NamedFigure(TitleStr);
  set(h, 'WindowStyle', 'docked');
  hold off
  plot(f, C, 'b.')
  if Freq > 0
    hold on
    CRange = [min(C), MaxC];
    plot([f(IndStart), f(IndStart)], CRange, 'g-')
    plot([f(IndStop), f(IndStop)], CRange, 'g-')
    plot(Freq, MaxC, 'ro', 'MarkerFaceColor', 'r')
    hold off
    fStop = max([2 * Freq, f(IndStart) * 1.1]);
    xlim([0, fStop])
  else
    xlim([0, 5])
  end
  xlabel('Frequency (Hz)', 'FontSize', 18);
  ylabel('AutoCorrelation', 'FontSize', 18);
  title(RealUnderscores(TitleStr), 'FontSize', 18);
end

if MaxC < CCutoff
  Freq = 0;
end

%Calculate NumSigma if needed
if nargout > 1
  if MaxC <= CCutoff
    NumSigma = 0;
  else
    NumSigma = abs(MaxC - mean(C)) / std(C);
  end
  if nargout == 2
    varargout = {NumSigma};
  else
    varargout = {NumSigma, MaxC};
  end
else
  varargout = {};
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotVar = DoPlot(PlotSubject)
if ischar(PlotSubject)
  PlotVar = true;
else
  PlotVar = PlotSubject;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PeakInd, PeakVal, LookStart, LookStop] = FindSlowWavePeak(C)
DerivLen = 10;

PeakInd = [];
PeakVal = [];
LookStart = [];
LookStop = [];

StartInd = find(C <= 0, 1);  %insist on passing apoint of C = 0
if isempty(StartInd)
  return
end

Deriv = C((DerivLen + StartInd):end) - C(StartInd:(end-DerivLen));
Trivial = 0.2 * median(abs(Deriv));
Neg1 = find(Deriv < -Trivial, 1);
if isempty(Neg1)
  return
end
Pos1 = find(Deriv(Neg1:end) > Trivial, 1) + Neg1 - 1;
if isempty(Pos1)
  return
end
Neg2 = find(Deriv(Pos1:end) < -Trivial, 1) + Pos1 - 1;
if isempty(Neg2)
  return
end
Pos2 = find(Deriv(Neg2:end) > Trivial, 1) + Neg2 - 1;
if isempty(Pos2)
  Pos2 = length(Deriv);
end

LookStart = Pos1 + StartInd - 1;
LookStop = Pos2 + StartInd - 1;
%A peak should exist near Neg2 + StartInd - 1.
TestPeakInd = StartInd + Neg2 - 1;
LastNeg = find(C(TestPeakInd:end) <= 0, 1);  %Try to find another C <= 0
if ~isempty(LastNeg)
  LastNeg = LastNeg + TestPeakInd - 1;
  if(LastNeg > LookStop)
    LookStop = LastNeg;
  end
end
						     
[PeakVal, PeakInd] = max(C(LookStart:LookStop));
PeakInd = PeakInd + LookStart - 1;
return