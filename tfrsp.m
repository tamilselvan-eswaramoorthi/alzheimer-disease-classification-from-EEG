function [tfr,stft,time,f] = tfrsp(signal,time,len,win,trace);


if (nargin == 0),
 error('At least 1 parameter required');
end;
[xrow,xcol] = size(signal);
if (nargin < 1),
 error('At least 1 parameter is required');
elseif (nargin <= 2),
 len=xrow;
end;

hlength=floor(len/4);
hlength=hlength+1-rem(hlength,2);

if (nargin == 1),
 time=1:xrow; 
 win = tftb_window(hlength); 
 trace=0;
elseif (nargin == 2) | (nargin == 3),
 win = tftb_window(hlength); trace=0;
elseif (nargin == 4),
 trace = 0;
end;

if (len<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(time);
if (xcol==0) | (xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(len)~=len),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

[hrow,hcol]=size(win); Lh=(hrow-1)/2;
if (hcol~=1)|(rem(hrow,2)==0),
 error('H must be a smoothing window with odd length');
end;

tfr= zeros (len,tcol) ;  
if trace, disp('Spectrogram'); end;
for icol=1:tcol,
 ti= time(icol); tau=-min([round(len/2)-1,Lh,ti-1]):min([round(len/2)-1,Lh,xrow-ti]);
 indices= rem(len+tau,len)+1; 
 if trace, disprog(icol,tcol,10); end;
 tfr(indices,icol)=signal(ti+tau).*conj(win(Lh+1+tau))/norm(win(Lh+1+tau));
end;
stft=fft(tfr);
tfr=abs(fft(tfr)).^2; 
if trace, fprintf('\n'); end;

if (nargout==0),
 tfrqview(tfr,signal,time,'tfrsp',win);
elseif (nargout==4),
 if rem(len,2)==0, 
  f=[0:len/2-1 -len/2:-1]'/len;
 else
  f=[0:(len-1)/2 -(len-1)/2:-1]'/len;  
 end;
end;

