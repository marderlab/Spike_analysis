function abfStruct = LoadAbf(fileName, varargin)
% abfStruct = LoadAbf(fileName, wantedString1, wantedString2, ...)
% This script is a wrapper for abfload.m
% It organizes the data in an .abf file into a structure with
% fields:
%  header: structure with header information, as from load_abf
%  time: array of times (in ms)
%  data: structure with fields that are arrays, traces (voltage or
%        current) as stored in the .abf file
%  units: structure with fields that are strings, specifying the
%        units in data.  data and units have the same field names.

[data, dt_us, header]=abfload(fileName);
% Check to make sure data is numSamples x numChannels x numEpisodes
numSamples = header.dataPtsPerChan;
numChannels = header.nADCNumChannels;

switch header.nOperationMode
  case 3,
    if isfield(header, 'lEpisodesPerRun')
      numEpisodes = header.lEpisodesPerRun;
    elseif isfield(header, 'lActualEpisodes')
      numEpisodes = header.lActualEpisodes;
    else
      warning('WAVEFORM:LoadAbf', ...
             'Unable to determine numEpisodes from header, assuming only one.')
      numEpisodes = 1
    end
  case 5, numEpisodes = header.lActualEpisodes;
  otherwise, error('Currently only supports gap-free or episodic data.')
end
if numEpisodes == 0
  % patch for ATF->ABF files
  numEpisodes = 1;
end
[d1, d2, d3] = size(data);
if any([numSamples, numChannels, numEpisodes] ~= [d1, d2, d3])
  error('Data size doesn''t match header reading %s', fileName)
end

% save the sample times and header into the structure
sampleTimes = (1.0e-3 * dt_us) * (0:(numSamples-1)); % convert to ms
abfStruct.header = header;
abfStruct.time = sampleTimes;

% save the data and units of the desired fields
numFound = 0;
for m = 1:numChannels
  channelName = header.recChNames{m};
  channelUnits = header.recChUnits{m};
  
  if nargin > 1
    notWanted = true;
    for n = 1:length(varargin)
      if ~isempty(strfind(channelName, varargin{n}))
        notWanted = false;
        break
      end
    end
    if notWanted
      continue
    end
  end
  
  numFound = numFound + 1;
  channelData = squeeze(data(:,m,:));  
  channelName = replaceBadChars(channelName);
  abfStruct.data.(channelName) = channelData;
  abfStruct.units.(channelName) = channelUnits;
end

if numFound == 0 || (nargin > 0 && numFound < nargin)
  error('Unable to find requested wave data')
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fieldName = replaceBadChars(fieldName)
badCharInds = regexp(fieldName, '\W');
if ~isempty(badCharInds)
  fieldName(badCharInds) = '_';
end
return
