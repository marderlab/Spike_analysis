function varargout = Toc
s = toc;
if(nargout == 1)
  varargout = {s};  
else
  varargout = {};
  if(s < 60)
    TimeString = sprintf('%.6fs', s);
    disp(['Elapsed time: ', TimeString])
    return
  end
  m = floor(s/60);
  s = s - m * 60;
  if(m < 60)
    TimeString = sprintf('%gm%.6fs', m, s);
    disp(['Elapsed time: ', TimeString])
    return
  end
  h = floor(m/60);
  m = m - h * 60;
  if(h < 24)
    TimeString = sprintf('%gh%gm%.6fs', h, m, s);
    disp(['Elapsed time: ', TimeString])
    return
  end
  d = floor(h/24);
  h = h - d * 24;
  TimeString = sprintf('%gd%gh%gm%.6fs', d, h, m, s);
  disp(['Elapsed time: ', TimeString])
end
return