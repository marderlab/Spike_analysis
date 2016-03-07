function TestWaveform(whichTest, varargin)

if nargin < 1
  whichTest = 'all';
end

saveResults = any(strcmp(varargin, 'save'));
forceDisplay = any(strcmp(varargin, 'display'));

thisFile = mfilename('fullpath');
dataDir = [thisFile(1:(end-12)), 'Datasets/'];




if strcmp(whichTest, 'Jon') || strcmp(whichTest, 'all')
  warning('WAVEFORM:compareAnalyzeWaveform', ...
          'Implement AnalyzeWaveform specific tests')
  testWaveformDirectory('Jon', dataDir, saveResults, forceDisplay)
end

if strcmp(whichTest, 'Tilman') || strcmp(whichTest, 'all')
  warning('WAVEFORM:compareAnalyzeWaveform', ...
          'Implement AnalyzeWaveform specific tests')
  testWaveformDirectory('Tilman', dataDir, saveResults, forceDisplay)
end

if strcmp(whichTest, 'Maria') || strcmp(whichTest, 'all')
  warning('WAVEFORM:compareGetFICurve', 'Implement F-I specific tests')
  testWaveformDirectory('Maria', dataDir, saveResults, forceDisplay)
end

if strcmp(whichTest, 'Lamont') || strcmp(whichTest, 'all')
  testWaveform_Lamont(dataDir)
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dirName, baseName, fileType] = parseFileName(fileName)
ind = strfind(fileName, '/');
if any(ind)
  ind = ind(end);
  baseName = fileName((ind+1):end);
  dirName = fileName(1:(ind-1));
else
  baseName = fileName;
  dirName = '';
end

ind = strfind(baseName, '.');
if any(ind)
  ind = ind(end);
  fileType = baseName((ind+1):end);
  baseName = baseName(1:(ind-1));
else
  fileType = '';
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataList = getDataFiles(testDir)
% return list of all the data files in the given test directory
dirList = dir(testDir);
dataList = {};
for n = 1:length(dirList)
  name_n = dirList(n).name;
  [~, ~, fileType] = parseFileName(name_n);
  if strcmp(fileType, 'abf') || strcmp(fileType, 'mat')
    dataList = [dataList, [testDir, '/', name_n]]; %#ok<AGROW>
  end
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadData(dataFile)
[~, ~, fileType] = parseFileName(dataFile);
if strcmp(fileType, 'abf')
  assignin('caller', 'abfStruct', LoadAbf(dataFile))
else
  % assume matlab file format
  matData = load(dataFile, '-mat');
  varNames = fieldnames(matData);
  for n = 1:length(varNames)
    varName = varNames{n};
    assignin('caller', varName, matData.(varName))
  end
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [missing, multiple] = checkSpikeOverlap(spikesA, spikesB)

n1A = spikesA.n1List;
n2A = spikesA.n2List;
tA = spikesA.times;
n1B = spikesB.n1List;
n2B = spikesB.n2List;
tB = spikesB.times;
missing = [];
multiple = [];
for n = 1:length(n1A)
  % find spikes that overlap with the nth A-spike
  m = find(n1B < n2A(n) & n2B > n1A(n));
  % exclude spikes that are closer to (n-1)th or (n+1)th A-spike
  exclude = [];
  for k = 1:length(m)
    m_k = m(k);
    if (n > 1 && tB(m_k) - tA(n-1) < tA(n) - tB(m_k)) || ...
       (n < length(n1A) && tA(n+1) - tB(m_k) < tB(m_k) - tA(n))
      exclude = [exclude, m_k]; %#ok<AGROW>
    end
  end
  if any(exclude)
    m = setdiff(m, exclude);
  end

  % exclude spikes that are closer to (n+1)th spike
  
  
  if isempty(m)
    missing = [missing, n]; %#ok<AGROW>
  elseif length(m) >= 2
    
    multiple = [multiple, n]; %#ok<AGROW>
  end
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compareOK = compare(testCommand, results, resultsCorrect)
[compareCommand, ~] = strtok(testCommand, '(');
compareCommand = eval(['@compare', compareCommand]);
compareOK = compareCommand(results, resultsCorrect);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compareOK = compareGetSpikes(results, resultsCorrect)

[missingSpikes, multipleA] = checkSpikeOverlap(resultsCorrect, results);
if any(missingSpikes)
  missingTimes = resultsCorrect.times(missingSpikes) / 1000;
  disp('Missing spike(s) at time:')
  disp(missingTimes)
end
if any(multipleA)
  multipleTimes = resultsCorrect.times(multipleA) / 1000;
  disp('Multiple spikes found when there should be only one at time:')
  disp(multipleTimes)
end
[extraSpikes, multipleB] = checkSpikeOverlap(results, resultsCorrect);
if any(extraSpikes)
  extraTimes = results.times(extraSpikes) / 1000;
  disp('Extra spike(s) at time:')
  disp(extraTimes)
end
if any(multipleB)
  multipleTimes = results.times(multipleB) / 1000;
  disp('One spike found when there should be many at time:')
  disp(multipleTimes)
end

compareOK = ~(any(missingSpikes) ||...
              any(multipleA) || ...
              any(extraSpikes) || ...
              any(multipleB));
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compareOK = compareAnalyzeWaveform(results, resultsCorrect) %#ok<DEFNU>
compareOK = compareGetSpikes(results.spikes, resultsCorrect.spikes);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compareOK = compareGetFICurve(results, resultsCorrect) %#ok<DEFNU>
compareOK = compareGetSpikes(results.spikes, resultsCorrect.spikes);
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayNotes()
try
  notes = evalin('caller', 'notes');
  fprintf(2, ' ->notes: ');
  disp(notes)  
catch err
  if ~strcmp(err.identifier, 'MATLAB:UndefinedFunction')
    rethrow err
  end
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testWaveformDirectory(subDirName, dataDir, ...
                               saveResults, forceDisplay)
% test analysis of electrophysiology data in a given subdirectory

fprintf('Testing %s\n', subDirName)

testDir = [dataDir, subDirName];
dataFiles = getDataFiles(testDir);

for n = 1:length(dataFiles)
  name_n = dataFiles{n};
  [~, baseName, ~] = parseFileName(name_n);
  
  clear notes
  loadData(name_n)
  fprintf(2, '%s: ', baseName);
  
  saveName = ([testDir, '/', baseName, '.save_test']);
  if saveResults
    if forceDisplay
      resultsCorrect = eval(dispCommand);
    else
      resultsCorrect = eval(testCommand);
    end
    save(saveName, 'resultsCorrect')
    disp('Results saved.')
  elseif exist(saveName, 'file')
    loadData(saveName)
    if forceDisplay
      results = eval(dispCommand);
    else
      results = eval(testCommand);
    end
    compareOK = compare(testCommand, results, resultsCorrect);
    if compareOK
      disp('Passed.')
      displayNotes()
    else
      disp('Saved results do not match, graphically displaying result.')
      if ~forceDisplay
        eval(dispCommand)
      end
      displayNotes()
      keyboard
    end
  else
    disp('Saved results do not exist, graphically displaying result.')
    eval([dispCommand, ';']);
    displayNotes()
  end
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testWaveform_Lamont(dataDir)

testDir = [dataDir, 'Lamont'];
extracellularFile = [testDir, '/', '753_079_0036_pdn_extra.abf'];

[t, pdn] = GetExtracellular(extracellularFile, 'pdn');

extra = AnalyzeExtracellular(t, pdn(1,:), 'plotSubject', 'pdn', ...
                             'debugPlots', true);
extra.burst
extra.burst.durations

return