function [QTFC,C,wind,AF,paramOUT,LOG] = TFCspwv(x1,x2,sig_len,vec,tau,lambda,fs);


x1 = hilbert(real(x1));
x2 = hilbert(real(x2));
% auto spectra
[C(:,:,1),AF(:,:,1),wind] = TFspwv(x1(:),sig_len,vec,tau,lambda);
[C(:,:,2),AF(:,:,2)] = TFspwv(x2(:),sig_len,vec,tau,lambda);
% cross spectra
wvx = tfrwv([x1(:),x2(:)] ,1:length(x1),sig_len);
A12=fft(ifft(wvx).');clear wvx
Asmooth12=( A12.*fftshift(wind) ).';
AF(:,:,3) = A12;clear A12
C(:,:,3)= ifft(fft(Asmooth12).').';clear Asmooth12
clear A12
% Quadratic coherence
QTFC= ( (C(:,:,3).*conj(C(:,:,3))) ./ (C(:,:,1).*C(:,:,2)) );

if nargout>4
    log(1) = sum(sum(QTFC>1));
    log(2) = sum(sum(QTFC<0));
    log2(1) = log(1)/(size(QTFC,1)*size(QTFC,2)) ;
    log2(2) = log(2)/(size(QTFC,1)*size(QTFC,2)) ;
    [fwhm_t,fwhm_f,Aperc_t,Aperc_f,wind_t,wind_f] = TFspwv_res(x1(:),sig_len,vec,tau,lambda,90,fs);
    paramOUT.out = log;
    paramOUT.out_perc = log2;
    paramOUT.fwhm_t = fwhm_t;
    paramOUT.fwhm_f = fwhm_f;
    paramOUT.Area90perc_t = Aperc_t;
    paramOUT.Area90perc_f = Aperc_f;
    paramOUT.wind_t = wind_t;
    paramOUT.wind_f = wind_f;
    paramOUT.t = [1:length(x1)]/fs;
    paramOUT.F = linspace(0,fs/2,sig_len);
    paramOUT.v0 = vec;
    paramOUT.tau0 = tau;
    paramOUT.lambda = lambda;
    paramOUT.Fbins = sig_len;
    
end
if nargout>5
    log_leg=char('QTFC>1 ','QTFC<0 ','Dt(fwhm)','Dt(90%)','Df(fwhm)','Df(90%)');
    LOG = [(log_leg),[' ';' ';' ';' ';' ';' '],num2str([log(:);[round(fwhm_t*fs);round(Aperc_t*fs);round(fwhm_f*sig_len/(fs/2));round(Aperc_f*sig_len/(fs/2))]]),[' samp  ';' samp  ';' samp  ';' samp  ';' samp  ';' samp  '] ...
        num2str([log2(:);[fwhm_t;(Aperc_t);fwhm_f;Aperc_f]],3),[' %  ';' %  ';' sec';' sec';' Hz ';' Hz '] ];
end