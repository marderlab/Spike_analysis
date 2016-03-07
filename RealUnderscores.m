function TexSafeString = RealUnderscores(TitleStr)
% TexSafeString = RealUnderscores(TitleStr)
%  Protects the underscores in TitleStr from being interpreted as
%  subscript commands by inserting a '\' before them.

Ind = strfind(TitleStr, '_');
TexSafeString = '';
Previous = 1;
for n = 1:length(Ind)
  if(Ind(n) > 1)
    TexSafeString = [TexSafeString, TitleStr(Previous:(Ind(n)-1))];
  end
  TexSafeString = [TexSafeString, '\'];
  Previous = Ind(n);
end
TexSafeString = [TexSafeString, TitleStr(Previous:end)];
return