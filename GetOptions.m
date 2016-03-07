function [options, varargout] = GetOptions(defaults, args, allowExtra)
% [options, modified] = GetOptions(defaults, args)
% Create a structure of options from defaults and user-overrides.
% INPUTS:
%  -defaults: a cell array with keyword, value ordering
%             or a structure
%  -args: either a structure or
%                a cell array
%         if a cell array, user can pass ALL keywords or
%                                   just go by input order.
%         if args is ommitted or length zero, the options will be the defaults.
%  OPTIONAL
%  -allowExtra    false
%    if set to true, ignore options specified in args that don't correspond
%    to anything in defaults.
% OUTPUTS:
%  -options: structure with fields set to their proper values

if nargin < 1 || nargin > 3
  help GetOptions
  error('Invalid number of inputs.')
end
if ~iscell(defaults) && ~isstruct(defaults)
  help GetOptions
  error('"defaults" must be a cell array or structure.')
end
if iscell(defaults);
  try
    defaults = cellOptionsToStruct(defaults);
  catch err
    % get help for the m-file that called GetOptions
    getHelp()
    % rethrow error
    rethrow(err)
  end
end


fNames = fieldnames(defaults);
modified = struct();
for n = 1:length(fNames)
  modified.(fNames{n}) = false;
end

if nargin == 1 || isempty(args)
  options = defaults;
  varargout = {modified};
  return
end

if nargin < 3
  allowExtra = false;
end

if iscell(args) && length(args) == 1 && isstruct(args{1})
  args = args{1};
end

if isstruct(args)
  % user passed a struct containing options
  try
    [options, modified] = ...
      getStructOptions(defaults, modified, args, allowExtra);
  catch err
    % get help for the m-file that called GetOptions
    getHelp()
    % rethrow error
    rethrow(err)
  end
elseif iscell(args)
  % user passed a cell array containing options
  try
    [options, modified] = ...
      getCellArrOptions(defaults, modified, args, allowExtra);
  catch err
    % get help for the m-file that called GetOptions
    getHelp()
    % rethrow error
    rethrow(err)
  end
else
  help GetOptions
  error('"args" must be a structure or a cell array.')
end
varargout = {modified};
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optStruct = cellOptionsToStruct(cellOpts)
numOpts = length(cellOpts);
if mod(numOpts, 2) ~= 0
  error('options must come in "key" "value" pairs.')
end
numOpts = numOpts / 2;
for n = 1:numOpts
  mVal = 2 * n;
  mKey = mVal - 1;
  key = cellOpts{mKey};
  val = cellOpts{mVal};
  if ~ischar(key)
    error('Invalid option key: %s', key)
  end
  optStruct.(key) = val;
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [options, modified] = ...
  getStructOptions(defaults, modified, argStruct, allowExtra)
% get options as defaults with overrides specified by argStruct

options = defaults;
fNames = fieldnames(argStruct);
for n = 1:length(fNames)
  field_n = fNames{n};
  value_n = argStruct.(field_n);
  if isempty(value_n)
    % skip empty values
    continue
  end
  
  if isfield(options, field_n) || allowExtra
    % this is a valid field, override the default
    %              or
    % save extra field (maybe another program will want it)
    options.(field_n) = value_n;
    modified.(field_n) = true;
  else
    error('Invalid option: %s', field_n)
  end
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [options, modified] = ...
  getCellArrOptions(defaults, modified, args, allowExtra)
% get options as default with overrides specified by cell array args

fNames = fieldnames(defaults);
keys = args(1:2:end);
options = defaults;

if mod(length(args), 2) == 0 && iscellstr(keys) && ...
   (~isempty(intersect(fNames, keys)) || length(args) > length(fNames))
  % overrides are keyword, value pairs
  
  if ~allowExtra && length(keys) > length(fNames)
    error('Too many arguments')
  end
  values = args(2:2:end);
  for n = 1:length(keys)
    if isempty(values{n})
      % skip empty values
      continue
    end
    
    if isfield(defaults, keys{n}) || allowExtra
      % this is a valid field, override the default
      %              or
      % save extra field (maybe another program will want it)
      options.(keys{n}) = values{n};
      modified.(keys{n}) = true;
    else
      error('Invalid option: %s', keys{n});
    end
  end
else
  % overrides are just a list of values
  
  if length(args) > length(fNames)
    error('Too many arguments')
  end
  
  for n = 1:length(args)
    options.(fNames{n}) = args{n};
    modified.(fNames{n}) = true;
  end
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getHelp()
% get help for the m-file that called GetOptions

st = dbstack;
thisFile = mfilename;
for n = 1:length(st)
  if ~strcmp(st(n).file, thisFile)
    %this is the calling file
    callFile = st(n).file;
    ind = find(callFile == filesep, 1, 'last');
    if ~isempty(ind)
      callFile = callFile((ind+1):end);
    end
    help(callFile)
    return
  end
end

help GetOptions
return
