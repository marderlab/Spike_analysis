function varargout=f(dt,x,nw,K,W_keep,conf_level,p_FFT_extra,tapers)


% dt is a scalar

% nw is the desired time-bandwidth product.  The frequency resolution is
%   given by nw/(N*dt).  Usually, nw=4 is a good place to start.
% conf_level is the confidence level of computed confidence intervals
% N_fft is the length to which data is zero-padded before FFTing
%
% f is the frequency base, which is one-sided
% the varargouts are the sigmas

% this version uses fft(), but does it on one sample, and one taper, at a
% time.  Also, it only stores the output up to frequency W_keep.  These 
% changes make it much more space-efficient that the "standard" multitaper
% code that I wrote.  Also, it's very fast.

% this version calculates the autocorrelation, not the autocovariance.
% I.e. we don't subtract off the mean first

% this version works with an arbitrary number of signals

% get the timing info, calc various scalars of interest
N=size(x,1);  % number of time points per process sample
Q=size(x,2);  % number of signals
R=size(x,3);  % number of samples of the process
fs=1/dt;

% process args
if nargin<4 || isempty(K)
  K=2*nw-1;
end
if nargin<5 || isempty(W_keep)
  W_keep=fs/2;
end
if nargin<6 || isempty(conf_level)
  conf_level=0;
end
if nargin<7 || isempty(p_FFT_extra)
  p_FFT_extra=2;
end
if nargin<8
  tapers=[];
end

% N for FFT
N_fft=2^(ceil(log2(N))+p_FFT_extra);

% compute frequency resolution
f_res_diam=2*nw/(N*dt);

% generate the dpss tapers if necessary
persistent N_memoed nw_memoed K_memoed tapers_memoed;
if isempty(tapers) 
  if isempty(tapers_memoed) | ...
     N_memoed~=N | nw_memoed~=nw | K_memoed~=K
    %fprintf(1,'calcing dpss...\n');
    tapers_memoed=dpss(N,nw,K);
    N_memoed=N;
    nw_memoed=nw;
    K_memoed=K;
  end
  tapers=tapers_memoed;
end
tapers=reshape(tapers,[N 1 1 K]);

% generate the frequency base
% hpfi = 'highest positive frequency index'
hpfi=ceil(N_fft/2);
f=fs*(0:(hpfi-1))'/N_fft;
f=f(f<=W_keep);
N_f=length(f);

% taper and do the FFTs
if N_f*Q*R*K<1e5
  % if dimensions are not too big, do this the easy way
  x_tapered=repmat(tapers,[1 Q R 1]).*repmat(x,[1 1 1 K]);
  X=fft(x_tapered,N_fft);
  X=X(1:N_f,:,:,:);
else
  % if dimensions are big, do this in a more space-efficient way
  X=zeros([N_f Q R K]);
  for r=1:R  % windows
    for k=1:K  % tapers
      x_this_tapered=repmat(tapers(:,:,:,k),[1 Q 1]).*x(:,:,r);
      X_this=fft(x_this_tapered,N_fft);
      X(:,:,r,k)=X_this(1:N_f,:);
    end
  end
end

% % convert to power by squaring, and to a density by dividing by fs
% Pxxs=(abs(X).^2)/fs;
% Pyys=(abs(Y).^2)/fs;
% Pyxs=(Y.*conj(X))/fs;

% need to generate all the cross-power spectra, with auto-power spectra on
% the diagonal
Pxys=zeros([N_f Q Q R K]);
for qj=1:Q
  conjXqj=conj(X(:,qj,:,:));
  for qi=1:Q
    Pxys(:,qi,qj,:,:)=X(:,qi,:,:).*conjXqj;
  end
end
Pxys=Pxys/fs;  % convert to density

% multiply by 2 (i.e. make into one-sided power spectra)
Pxys=2*Pxys;

% plot stuff
%for j=1:K
%  n_plots=6;
%  figure;
%  subplot(n_plots,1,1);
%  plot(t,x_tapered(:,:,j));
%  ylabel('x tapered');
%  title(sprintf('Taper %d',j));
%  subplot(n_plots,1,2);
%  plot(t,y_tapered(:,:,j));
%  ylabel('y tapered');
%  subplot(n_plots,1,3);
%  plot(f,Pxxs(:,:,j));
%  ylabel('Pxx');
%  subplot(n_plots,1,4);
%  plot(f,Pyys(:,:,j));
%  ylabel('Pyy');
%  subplot(n_plots,1,5);
%  plot(f,abs(Pyxs(:,:,j)));
%  ylabel('abs Pyx');
%  subplot(n_plots,1,6);
%  plot(f,angle(Pyxs(:,:,j)));
%  ylabel('angle Pyx');
%end

% _sum_ across samples, tapers (keep these around in case we need to 
% calculate the take-away-one spectra for error bars)
PxyRK=sum(sum(Pxys,5),4);
% PxyRK is of shape [N_f Q Q]

% convert the sum across samples, tapers to an average; these are our 
% 'overall' spectral estimates
Pxy=PxyRK/(R*K);
% Pxy is of shape [N_f Q Q]

% plot stuff
%n_plots=6;
%figure;
%subplot(n_plots,1,1);
%plot(t,x);
%ylabel('x');
%title(sprintf('Final',j));
%subplot(n_plots,1,2);
%plot(t,y);
%ylabel('y');
%subplot(n_plots,1,3);
%plot(f,Pxx);
%ylabel('Pxx');
%subplot(n_plots,1,4);
%plot(f,Pyy);
%ylabel('Pyy');
%subplot(n_plots,1,5);
%plot(f,abs(Pyx));
%ylabel('abs Pyx');
%subplot(n_plots,1,6);
%plot(f,angle(Pyx));
%ylabel('angle Pyx');

% % calculate coherence
% Cyx=Pyx./sqrt(Pxx.*Pyy);

% calculate coherence
Pxx=zeros([N_f Q 1]);
for q=1:Q
  Pxx(:,q)=Pxy(:,q,q);
end
Pyy=reshape(Pxx,[N_f 1 Q]);
Cxy=Pxy./sqrt(repmat(Pxx,[1 1 Q]).*repmat(Pyy,[1 Q 1]));

% separate out magnitude, phase
Cxy_mag=abs(Cxy);
Cxy_phase=unwrap(angle(Cxy));

% calc the sigmas
if conf_level>0
  % calculate the transformed power, coherence magnitude
  Pxx_xf=log10(Pxx);
  Cxy_mag_xf=atanh(Cxy_mag);

  % calculate the take-away-one spectra
  Pxxs=zeros([N_f Q]);
  PxxRK=zeros([N_f Q]);
  for q=1:Q
    Pxxs(:,q)=Pxys(:,q,q);    
    PxxRK(:,q)=PxyRK(:,q,q);
  end 
  Pxys_tao=(repmat(PxyRK,[1 1 1 R K])-Pxys)/(R*K-1);
  
  % calc the take-away-one coherence
  Pxxs_tao=zeros([N_f Q 1 R K]);
  for q=1:Q
    Pxxs_tao(:,q,:,:,:)=Pxys_tao(:,q,q,:,:);
  end
  Pyys_tao=reshape(Pxxs_tao,[N_f 1 Q R K]);
  Cxys_tao=Pxys_tao./sqrt(repmat(Pxxs_tao,[1 1 Q 1 1]).*...
                          repmat(Pyys_tao,[1 Q 1 1 1]));

  % transform the take-away-one spectra, coherence
  Pxxs_tao_xf=log10(Pxxs_tao);
  Cxys_tao_mag=abs(Cxys_tao);
  Cxys_tao_mag_xf=atanh(Cxys_tao_mag);
  Cxys_tao_phase=angle(Cxys_tao);

  % calculate the sigmas on the spectra
  Pxxs_tao_xf_mean=mean(mean(Pxxs_tao_xf,5),4);
  Pxx_xf_sigma=...
    sqrt((R*K-1)/(R*K)*...
         sum(sum((Pxxs_tao_xf-...
                  repmat(Pxxs_tao_xf_mean,[1 1 1 R K])).^2,5),4));

  % calculate the coherence magnitude sigma
  Cxys_tao_mag_xf_mean=mean(mean(Cxys_tao_mag_xf,5),4);
  Cxy_mag_xf_sigma=...
    sqrt((R*K-1)/(R*K)*...
         sum(sum((Cxys_tao_mag_xf-...
                  repmat(Cxys_tao_mag_xf_mean,[1 1 1 R K])).^2,5),4));

  % calculate the coherence phase sigma
  Cxys_tao_hat=Cxys_tao./Cxys_tao_mag;
  Cxys_tao_hat_mean=mean(mean(Cxys_tao_hat,5),4);
  Cxy_phase_sigma=sqrt(2*(R*K-1)*(1-abs(Cxys_tao_hat_mean)));

  % calculate the confidence intervals
  ci_factor=tinv((1+conf_level)/2,R*K-1);
  Pxx_ci(:,:,1)=10.^(Pxx_xf-ci_factor*Pxx_xf_sigma);
  Pxx_ci(:,:,2)=10.^(Pxx_xf+ci_factor*Pxx_xf_sigma);
  Cxy_mag_ci(:,:,:,1)=tanh(Cxy_mag_xf-ci_factor*Cxy_mag_xf_sigma);
  Cxy_mag_ci(:,:,:,2)=tanh(Cxy_mag_xf+ci_factor*Cxy_mag_xf_sigma);
  Cxy_phase_ci(:,:,:,1)=Cxy_phase-ci_factor*Cxy_phase_sigma;
  Cxy_phase_ci(:,:,:,2)=Cxy_phase+ci_factor*Cxy_phase_sigma;

  % assign the return values, returning sigmas
  varargout={f ...
             Pxx ...
             Cxy_mag Cxy_phase ...
             N_fft f_res_diam ...
             Pxx_ci ...
             Cxy_mag_ci ...
             Cxy_phase_ci ...
             Pxx_xf Pxx_xf_sigma ...
             Cxy_mag_xf Cxy_mag_xf_sigma ...
             Cxy_phase_sigma ...
             Pxxs_tao ...
             };
else
  % assign the return values, w/o CIs, sigmas
  varargout={f Pxx Cxy_mag Cxy_phase N_fft f_res_diam};
end
