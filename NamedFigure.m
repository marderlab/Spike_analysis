function varargout = NamedFigure(name)
% h = NamedFigure(name)
%  Looks for a figure with requested name, sets it to the current figure,
%  and returns a handle h to the figure.  If no such figure exists, it
%  creates one, names it, and returns a handle
%    INPUT:
%     -name   The name of the figure
%    OUTPUT:
%     -h      Handle to the figure

% find any figures with requested name
h = findobj('Name', name);

if isempty(h)
  % no such figure exists, create one
  
  % set some default properties for figures
  set(0, 'defaultaxesfontname', 'Arial')
  set(0, 'defaulttextfontname', 'Arial')
  set(0, 'defaultaxesfontsize', 15)
  set(0, 'defaulttextfontsize', 18)  
  
  % create the figure
  h = figure;
  
  % set convenient figure properties
  set(h, 'Name', name)
  set(h, 'RendererMode', 'manual')
  set(h, 'Renderer', 'painters')
  % Not much point to this, because it gets overwritten, but oh well:
  title(name);
elseif ishandle(h)
  % figure exists, get it
  set(0, 'CurrentFigure', h)
else
  % I've actually forgotten why this is here, presumably it causes no
  % trouble
  h = figure(h);
end

% return the figure handle if it's requested
switch nargout
 case 0, varargout = {};
 case 1, varargout = {h};
 otherwise, error('Too many output arguments.');
end
return