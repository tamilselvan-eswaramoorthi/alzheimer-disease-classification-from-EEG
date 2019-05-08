function [Q,C,paramOUT,LOG] = TFCmtsp(sig_one,sig_two,Fbins,wh,fs);
[hl,M]=size(wh);
if M>hl
    error('wh must be a column-vector')
end

for k = 1:M
    [spt,stft] = tfrsp(sig_one,[1:length(sig_one)],Fbins,wh(:,k)) ;
    STFT(:,:,1,k)=stft;
    
    [spt,stft] = tfrsp(sig_two,[1:length(sig_two)],Fbins,wh(:,k)) ;
    STFT(:,:,2,k)=stft;    
end

C=abs(STFT).^2;
C(:,:,3,:)=STFT(:,:,1,:).*conj(STFT(:,:,2,:));

for m=1:M
Sm1 = mean(C(:,:,1,1:m),4);
Sm2 = mean(C(:,:,2,1:m),4);
Sm12 = mean(C(:,:,3,1:m),4);
Q(:,:,m) = abs(Sm12).^2./(Sm1.*Sm2);
end

log(1) = sum(sum(Q(:,:,M)>1));
log(2) = sum(sum(Q(:,:,M)<0));
log2(1) = log(1)/(size(Q(:,:,M),1)*size(Q(:,:,M),2)) ;
log2(2) = log(2)/(size(Q(:,:,M),1)*size(Q(:,:,M),2)) ;
log_leg=char('QTFC>1 ','QTFC<0 ','Dt(fwhm)','Dt(90%)','Df(fwhm)','Df(90%)');
[fwhm_t,fwhm_f,Aperc_t,Aperc_f,wind_t,wind_f] = TFmtsp_res(wh,90,fs,0);

paramOUT.out = log;
paramOUT.out_perc = log2;
paramOUT.fwhm_t = fwhm_t;
paramOUT.fwhm_f = fwhm_f;
paramOUT.Area90perc_t = Aperc_t;
paramOUT.Area90perc_f = Aperc_f;
paramOUT.wind_t = wind_t;
paramOUT.wind_f = wind_f;
paramOUT.hf = wh;
paramOUT.t = [1:length(sig_one)]/fs;
paramOUT.F = linspace(0,fs,Fbins);
paramOUT.Fbins = Fbins;


LOG = [(log_leg),[' ';' ';' ';' ';' ';' '],num2str([log(:);[round(fwhm_t(end)*fs);round(Aperc_t(end)*fs);round(fwhm_f(end)*Fbins/fs);round(Aperc_f(end)*Fbins/fs)]]),[' samp  ';' samp  ';' samp  ';' samp  ';' samp  ';' samp  '] ...
num2str([log2(:);[fwhm_t(end);(Aperc_t(end));fwhm_f(end);(Aperc_f(end))]]),[' %  ';' %  ';' sec';' sec';' Hz ';' Hz '] ];
