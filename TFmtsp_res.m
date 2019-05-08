function [fwhm_t,fwhm_f,Aperc_t,Aperc_f,wind_t,wind_f] = TFmtsp_res(h,area_perc,fs,trace)



if nargin<4
    trace=0;
end
[hl,hr]=size(h);
if hr>hl
    error('windows must be a column vector or a matrix (length(h),numTaper)')
end

Fbins=2048;
t=[0:size(h,1)-1]/fs; F = [0:Fbins-1]/Fbins*fs;
for ih=1:size(h,2)
    wind = (h(:,1:ih)).^2;
    wind = mean(wind,2);
    wind=wind/max(wind);
    
    factsampT=0.1;
    samp2=linspace(t(1),t(end),length(t)/factsampT);
    wind2_int=interp1(t,wind,samp2);
    
    [d,ii2]=min(abs(samp2-t(end)/2)); % the max is not always the center
    AA(1)=wind2_int(ii2);
    jjend=min(ii2,length(samp2)-ii2)-1;
    for jj=2:jjend
        AA(jj)= AA(jj-1) + wind2_int(ii2-(jj-1)) + wind2_int(ii2+(jj-1));
    end
    ee=find(AA/sum(wind2_int)>area_perc/100);
    Aperc_t(ih)=ee(1)*2*abs(median(diff(samp2)));
    
    m05 = find(wind2_int>0.5);
    fwhm_t(ih)= abs(samp2(m05(1))-samp2(m05(end)));
    wind_t(:,ih)=fftshift(wind);
    clear wind samp*
    
    wind = fftshift(abs(fft(h(:,1:ih),Fbins)).^2);
    wind = (mean(wind,2));
    wind=wind/max(wind);
    
    factsampF=0.1;
    samp2=linspace(F(1),F(end),length(F)/factsampT);
    wind2_int=interp1(F,wind,samp2);
    
    [d,ii2]=min(abs(samp2-fs/2)); % the max is not always the center
    AA(1)=wind2_int(ii2);
    jjend=min(ii2,length(samp2)-ii2)-1;
    for jj=2:jjend
        AA(jj)= AA(jj-1) + wind2_int(ii2-(jj-1)) + wind2_int(ii2+(jj-1));
    end
    ee=find(AA/sum(wind2_int)>area_perc/100);
    Aperc_f(ih)=ee(1)*2*abs(median(diff(samp2)));
    
    m05 = find(wind2_int>0.5);
    fwhm_f(ih)= abs(samp2(m05(1))-samp2(m05(end)));
    wind_f(:,ih)=fftshift(wind);
    clear wind samp*
end

if trace==1
    display('colums represents the resolution of the average window')
    log_leg=char('Dt(fwhm)',['Dt(',num2str(area_perc),'%)'],'Df(fwhm)',['Df(',num2str(area_perc),'%)']);
    LOG = [(log_leg),[' ';' ';' ';' '],num2str([fwhm_t;Aperc_t;fwhm_f;Aperc_f]),[' sec';' sec';' Hz ';' Hz ']];
    display(LOG)
end

