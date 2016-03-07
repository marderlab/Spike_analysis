function Burst = AnalyzeBurst(Spike, SlowWave, t)
%   -Burst:  structure with burst information
%      -Burst.Freq is mean burst frequency (Hz)
%      -Burst.SpikeFreq is mean within-burst spike frequency (Hz)
%      -Burst.DutyCycle is the average burst duration/period
%      -Burst.Times is a plain list of burst times (ms)
%      -Burst.Durations is a list of burst durations (ms)
%      -Burst.SpikesPerBurst is a list of spikes per burst
%      -Burst.SpikeFrequencies is a list of spike frequencies (Hz)
%      -Burst.InterBurstIntervals is a list of intervals between bursts (ms)
%      -Burst.StartInds is a list of indices where bursts begin
%      -Burst.StopInds is a list of indices where bursts end

if Spike.freq <= 0 || (SlowWave.NumOnlySpike >= SlowWave.NumBurst)
  %This is not a bursting waveform
  Burst.Times = [];
  Burst.Durations = [];
  Burst.SpikesPerBurst = [];
  Burst.SpikeFrequencies = [];
  Burst.InterBurstIntervals = [];
  Burst.StartInds = [];
  Burst.StopInds = [];
  Burst.Freq = 0;
  Burst.SpikeFreq = 0;
  Burst.DutyCycle = 0;
  return
end

BurstInds = SlowWave.MinInds;
NumBurst = length(BurstInds) - 1;
CosPhase = cos(SlowWave.Phases);
t_Off = -Inf;

Burst.Times = zeros(NumBurst, 1);
Burst.Durations = zeros(NumBurst, 1);
Burst.SpikesPerBurst = zeros(NumBurst, 1);
Burst.SpikeFrequencies = zeros(NumBurst, 1);
Burst.InterBurstIntervals = zeros(NumBurst - 1, 1);
Burst.StartInds = zeros(NumBurst, 1);
Burst.StopInds = zeros(NumBurst, 1);

for BurstNum = 1:NumBurst
  Last_Off = t_Off;
  n1 = BurstInds(BurstNum);
  n2 = BurstInds(BurstNum + 1);
  %t1 = t(n1);
  %t2 = t(n2);
  n_On = find(CosPhase(n1:n2) < 0, 1) + n1 - 1;
  n_Off = find(CosPhase(n_On:n2) < 0, 1, 'last') + n_On - 1;
  t1 = t(n_On);
  t2 = t(n_Off);
  
  %SpikeInds = find(Spike.n1List > n1 & Spike.n2List <= n2);
  SpikeInds = find(Spike.times > t1 & Spike.times < t2);
  Spike_n1s = Spike.n1List(SpikeInds);
  Spike_n2s = Spike.n2List(SpikeInds);

  %{
  %This long block extends bursts to include all spikes that belong
  %in them
  if(length(SpikeInds) > 1 && ...
     (Spike_n1s(1) < n_On || Spike_n2s(end) > n_Off))
    %either expand burst to include spikes, or
    %  remove non-burst spikes from list
    InBurst = find(Spike_n2s >= n_On & Spike_n1s <= n_Off);
    if(length(InBurst) > 1)
      if(Spike_n1s(InBurst(1)) < n_On)
	n_On = Spike_n1s(InBurst(1));
      end
      if(Spike_n2s(InBurst(end)) > n_Off)
	n_Off = Spike_n2s(InBurst(end));
      end
      MaxDiff = max(diff(Spike_n1s(InBurst)));
      while(InBurst(1) > 1)
	NewDiff = Spike_n1s(InBurst(1)) - Spike_n2s(InBurst(1) - 1);
	if(NewDiff < MaxDiff)
	  InBurst = [InBurst(1)-1; InBurst(:)];
	  n_On = Spike_n1s(InBurst(1));
	  NewDiff = Spike_n1s(InBurst(2)) - n_On;
	  if(NewDiff > MaxDiff)
	    MaxDiff = NewDiff;
	  end
	else
	  break
	end
      end
      NumOrig = length(SpikeInds);
      while(InBurst(end) < NumOrig)
	NewDiff = Spike_n1s(InBurst(end)+1) - Spike_n2s(InBurst(end));
	if(NewDiff < MaxDiff)
	  InBurst = [InBurst(:); InBurst(end)+1];
	  n_Off = Spike_n2s(InBurst(end));
	  NewDiff = Spike_n1s(InBurst(end)) - Spike_n1s(InBurst(end-1));
	  if(NewDiff > MaxDiff)
	    MaxDiff = NewDiff;
	  end
	else
	  break
	end
      end      
    end
    SpikeInds = SpikeInds(InBurst);
  end
  %}
  t_On = t(n_On);
  t_Off = t(n_Off);
  NumSpikes = length(SpikeInds);
  Burst.StartInds(BurstNum) = n_On;
  Burst.StopInds(BurstNum) = n_Off;
  Burst.Times(BurstNum) = t_On;
  Burst.Durations(BurstNum) = t_Off - t_On;

  Burst.SpikesPerBurst(BurstNum) = NumSpikes;
  if(NumSpikes > 2)
    Interval = .001 * (Spike.times(SpikeInds(end)) ...
		       - Spike.times(SpikeInds(1)));  %sec
    Burst.SpikeFrequencies(BurstNum) = (NumSpikes-1) / Interval;
  else
    Burst.SpikeFrequencies(BurstNum) = NaN;
  end
  if(BurstNum > 1)
    Burst.InterBurstIntervals(BurstNum) = Last_Off - t_On;
  end
end

%If two bursts separate roughly continuous spiking, join them
%fprintf('NumBurst before refine: %g\n', NumBurst)
[Burst, NumBurst] = RefineBurst(Burst, Spike);
%fprintf('NumBurst after refine: %g\n', NumBurst)

if(NumBurst >= 3)
  Periods = .001 * (Burst.Times(2:end) - Burst.Times(1:(end-1)));
  Burst.Freq = mean(1.0 ./ Periods);
else
  Burst.Freq = 0;
end
Ind = find(isfinite(Burst.SpikeFrequencies));
Burst.SpikeFreq = mean(Burst.SpikeFrequencies(Ind));
Burst.DutyCycle = sum(Burst.Durations) / (t(BurstInds(end)) - t(BurstInds(1)));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Burst, NumBurst] = RefineBurst(Burst, Spike)
NumBurst = length(Burst.Times);
m = 1;
while(m < NumBurst)
  n1_m = Burst.StartInds(m);
  n2_m = Burst.StopInds(m);
  n1_mp1 = Burst.StartInds(m+1);
  n2_mp1 = Burst.StopInds(m+1);  
  
  mSpikeInds = find(Spike.n1List >= n1_m & Spike.n1List <= n2_m);
  mp1SpikeInds = find(Spike.n1List >= n1_mp1 & Spike.n1List <= n2_mp1);  
  mNum = length(mSpikeInds);
  mp1Num = length(mp1SpikeInds);
  %fprintf('\tmNum: %g, mp1Num: %g, DeltaInd: %g\n', ...
%	  mNum, mp1Num, mp1SpikeInds(1) - 1 - mSpikeInds(end))
  JoinBurst = false;
  if(mNum > 1 && mp1Num > 1 && mSpikeInds(end) + 1 == mp1SpikeInds(1))
    
    Spacing = Spike.times(mp1SpikeInds(1)) - Spike.times(mSpikeInds(end));
    mSTimes = Spike.times(mSpikeInds);
    mp1STimes = Spike.times(mp1SpikeInds);
    mSIntervals = diff(mSTimes);
    mp1SIntervals = diff(mp1STimes);
    MeanInt = 0.5 * (mean(mSIntervals) + mean(mp1SIntervals));
    STDInt = 0.5 * (std(mSIntervals) + std(mp1SIntervals));
    %fprintf('\t\tSpace: %g, MeanSpace: %g, 2 * STDSpace: %g\n', ...
    %        Spacing, MeanInt, 2 * STDInt)
    if(abs(MeanInt - Spacing) < 2 * STDInt)
      JoinBurst = true;
    end
  end
  
  if(JoinBurst)
    %fprintf('\t\t\tJoin!\n')
    Burst.Durations = [Burst.Durations(1:(m-1)); ...
	       Burst.Durations(m+1) + Burst.Times(m+1) - Burst.Times(m); ...
		       Burst.Durations((m+2):end)];
    Burst.SpikesPerBurst = [Burst.SpikesPerBurst(1:(m-1)); ...
		    Burst.SpikesPerBurst(m) + Burst.SpikesPerBurst(m + 1); ...
		    Burst.SpikesPerBurst((m+2):end)];
    DeltaT = mp1STimes(end) - mSTimes(1);
    Burst.SpikeFrequencies = [Burst.SpikeFrequencies(1:(m-1)); ...
		    (Burst.SpikesPerBurst(m)-1) / DeltaT; ...
		    Burst.SpikeFrequencies((m+2):end)];
    Burst.Times = [Burst.Times(1:m); Burst.Times((m+2):end)];
    Burst.StartInds = [Burst.StartInds(1:m); Burst.StartInds((m+2):end)];
    Burst.StopInds = [Burst.StopInds(1:(m-1)); Burst.StopInds((m+1):end)];
    Burst.InterBurstIntervals = [Burst.InterBurstIntervals(1:(m-1)); ...
		    Burst.InterBurstIntervals((m+1):end)];   
    NumBurst = NumBurst - 1;
  else
    m = m + 1;
  end
end

return
