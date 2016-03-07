function [t, varargout] = GetExtracellular(fileName, varargin)
% [t, ex1, ex2, ...] = GetExtracellular(fileName, ex1Name, ex2Name, ...)

abfS = LoadAbf(fileName);

t = abfS.time'; % (in ms)

fNames = fieldnames(abfS.units);
if isempty(varargin)
  % no names supplied, so return all extracellular data
  varargout = cell(1, length(fNames));
  for n = 1:length(fNames)
    varargout{n} = abfS.data.(fNames{n});
  end
else
  % get matching names
  numRequest = length(varargin);
  varargout = cell(1, numRequest);
  for n = 1:numRequest
    m = getMatchInd(varargin{n}, fNames);
    varargout{n} = abfS.data.(fNames{m})';
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matchInd = getMatchInd(requestName, fNames)
matchInd = [];
for n = 1:length(fNames)
  if StringCheck(fNames{n}, requestName)
    matchInd = [matchInd, n];
  end
end

if isempty(matchInd)
  fprintf(2, 'Available extracellulars:')
  for n = 1:length(fNames)
    fprintf(2, ' %s', fNames{n})
  end
  fprintf(2, '\n')
  error('Couldn''t find match for requested trace: %s', requestName)
elseif length(matchInd) > 1
  fprintf(2, 'Matching extracellulars:')
  for n = matchInd
    fprintf(2, ' %s', fNames{n})
  end
  fprintf(2, '\n')
  error('Found multiple matches requested trace: %s', requestName)  
end
return