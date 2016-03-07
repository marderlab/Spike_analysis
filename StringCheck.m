function isMatch = StringCheck(checkString, pattern, strictCase)
% isMatch = StringCheck(checkString, pattern, strictCase)
% Checks to see if checkString contains pattern.
%  INPUT PARAMETERS:
%   -checkString: string to be checked
%   -pattern: pattern to test against
%    OPTIONAL:
%    -strictCase:  If set to true, upper/lower case must be exact.
%           (Defaults to false)
%  OUTPUT PARAMETERS:
%   -isMatch:  boolean specifying if there is a match
if nargin < 3
  strictCase = false;
end

if ~strictCase
  checkString = lower(checkString);
  pattern = lower(pattern);
end

isMatch = ~isempty(strfind(checkString, pattern));
return