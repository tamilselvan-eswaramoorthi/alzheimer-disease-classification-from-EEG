function [fwhm_t,fwhm_f,Aperc_t,Aperc_f,wind_t,wind_f] = TFspwv_res(sig,rati_len,vec,tau,lambda,area_perc,freq,trace);



if nargin<8
    trace=0;
elseif nargin<7
    error('not sufficient inputs');
end

if area_perc<1
    error('area_perc ranges in [1,100]')
end

t=[0:length(sig)-1]/freq;
F=[0:rati_len-1]/rati_len*(freq/2);

Crow = length(sig);
Ccol = rati_len;

dy=[-Ccol/2:Ccol/2-1]/(Ccol/2);
if Crow/2- fix(Crow/2)==0
    dx=[-Crow/2:Crow/2-1]/(Crow/2);
else
    dx=[-Crow/2-1/2:Crow/2-3/2]/(Crow/2);
end


for il=1:length(lambda)
    for iv=1:length(vec)
        
        mu_v=  ((dx/vec(iv)).^2)  ;
        wind_v=exp(-pi*( mu_v.^2).^lambda(il) ); wind_v=wind_v/max(max(wind_v));
        wind=fftshift(abs(fft(wind_v)));
        wind=wind/max(wind);
        
        factsampT=.1;
        samp2=linspace(t(1),t(end),length(t)/factsampT);
        wind2_int=interp1(t,wind,samp2);
        
        [d,ii2]=max(wind2_int);
        AA(1)=wind2_int(ii2);
        
        jjend=min(ii2,length(samp2)-ii2)-1;
        for jj=2:jjend
            AA(jj)= AA(jj-1) + wind2_int(ii2-(jj-1)) + wind2_int(ii2+(jj-1));
        end
        ee=find(AA/sum(wind2_int)>area_perc/100);
        Aperc_t(il,iv)=ee(1)*2*abs(median(diff(samp2)));
        
        m05 = find(wind2_int>0.5);
        %         fwhm_t= abs(samp2(m05(1))-samp2(m05(end)));
        fwhm_t(il,iv)= abs(samp2(m05(1))-samp2(m05(end)));
        wind_t=wind(:);
        clear wind wind2 wind2_int
    end
    
    for itau=1:length(tau)
        mu_tau=  ((dy/tau(itau)).^2) ;
        wind_tau=exp(-pi*( mu_tau.^2).^lambda(il) );
        wind=fftshift(abs(fft(wind_tau)));
        wind=wind/max(wind);
        
        factsampF=.1;
        samp2=linspace(F(1),F(end),length(F)/factsampF);
        wind2_int=interp1(F,wind,samp2);
        
        [d,ii2]=max(wind2_int);
        AA(1)=wind2_int(ii2);
        
        jjend=min(ii2,length(samp2)-ii2)-1;
        for jj=2:jjend
            AA(jj)= AA(jj-1) + wind2_int(ii2-(jj-1)) + wind2_int(ii2+(jj-1));
        end
        ee=find(AA/sum(wind2_int)>area_perc/100);
        Aperc_f(il,itau)=ee(1)*2*abs(median(diff(samp2)));
        
        m05 = find(wind2_int>0.5);
        fwhm_f(il,itau)= abs(samp2(m05(1))-samp2(m05(end)));
        wind_f=wind(:);
    end
end

if trace
    display(['%------- FWHM-time [s] ------%'])
    display(num2str(fwhm_t,3))
    display(['%------- Aperc-time [s] ------%'])
    display(num2str(Aperc_t,3))
    
    display(['%------- FWHM-freq [Hz] ------%'])
    display(num2str(fwhm_f,3))
    display(['%------- Aperc-freq [Hz] ------%'])
    display(num2str(Aperc_f,3))
end
