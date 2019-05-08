function [tfr,timne,f] = tfrwv(sig,timne,len,trace);


if (nargin == 0),
    error('At least one parameter required');
end;
[xrow,xcol] = size(sig);

if (nargin == 1),
    timne=1:xrow; len=xrow ; trace=0;
elseif (nargin == 2),
    len=xrow ; trace=0;
elseif (nargin == 3),
    trace = 0;
end;

if (len<0),
    error('N must be greater than zero');
end;

[trow,tcol] = size(timne);
if (xcol==0)|(xcol>2),
    error('X must have one or two columns');
elseif (trow~=1),
    error('T must only have one row');
elseif (2^nextpow2(len)~=len),
    fprintf('For a faster computation, N should be a power of two\n');
end;

tfr= zeros (len,tcol);
if trace, disp('Wigner-Ville distribution'); end;
for icol=1:tcol,
    ti= timne(icol); taumax=min([ti-1,xrow-ti,round(len/2)-1]);
    tau=-taumax:taumax; indices= rem(len+tau,len)+1;
    tfr(indices,icol) = sig(ti+tau,1) .* conj(sig(ti-tau,xcol));
    tau=round(len/2);
    if (ti<=xrow-tau)&(ti>=tau+1),
        tfr(tau+1,icol) = 0.5 * (sig(ti+tau,1) * conj(sig(ti-tau,xcol))  + ...
            sig(ti-tau,1) * conj(sig(ti+tau,xcol))) ;
    end;
    if trace, disprog(icol,tcol,10); end;
end;
tfr= fft(tfr);
if (xcol==1), tfr=real(tfr); end ;

if (nargout==0),
    tfrqview(tfr,sig,timne,'tfrwv');
elseif (nargout==3),
    f=(0.5*(0:len-1)/len)';
end;
