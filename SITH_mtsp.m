function [sig_instre,tuning_param] = SITH_mtsp(varargin)


tuning_param.T = 816*4;    
tuning_param.Fbins = 2048; 
tuning_param.hfN = 901;    
tuning_param.hfM =4;        
tuning_param.hftm = 20;    
tuning_param.Ns = 100;     
tuning_param.fs = 4;       
tuning_param.statlev = [90,95,99];   
tuning_param.s_step = 5;   
tuning_param.doplot = 0;
pOUT_names = fieldnames(tuning_param);
J = 1;
while J <= length(varargin)
    key = varargin{J};
    keyN = find(strcmpi(key, pOUT_names));
    if isempty(keyN)
        warning('PPLIKEL:UnknownOption', 'Do not know option ''%s''', key);
    else
        tuning_param.(pOUT_names{keyN}) = varargin{J+1};
    end
    J = J + 2;
end
for id=1:length(pOUT_names)
     eval([char(pOUT_names(id)),'=paramOUT.',char(pOUT_names(id)),';'])
end

% Hermite functions
[hf] = hermf(hfN,hfM,hftm);
hf = hf.'; 
% Estimating TF resolution
[fwhm_t,fwhm_f] = TFmtsp_res(hf,90,fs);
display(['using hfN=',num2str(hfN),'; hftm=',num2str(hftm),'; hfM=',num2str(hfM),' tapers']);
display(['#surr =',num2str(Ns),' Fbins=',num2str(Fbins),' T=',num2str(T)]);
Dt = round(fwhm_t*fs);
Df = round(fwhm_f*Fbins/fs);
cTOT = nan(length(Df:Df:Fbins/4),length(Dt:Dt:T-Dt),Ns);
for is = 1:100

s1=RandStream('mrg32k3a','seed',is);
s2=RandStream('mrg32k3a','seed',is+Ns);
clear c2 c
% White noises
x1=hilbert(randn(s1,T,1));
x2=hilbert(randn(s2,T,1));
% Estimating SPWV-TFC
[c2]=TFCmtsp(x1,x2,Fbins,hf,fs);
% Reduce dimension
c = sqrt(c2(Df:Df:Fbins/4,Dt:Dt:end-Dt,end));
display(['iter=',num2str(is),' mean TFC=',num2str(mean(c(:)))]);
cTOT(:,:,is) =c;
end

for ii=1:round(size(cTOT,3)/s_step)
    for istat=1:length(statlev)
        SITH_TFC=prctile(cTOT(:,:,1:ii*s_step),statlev(istat),3);
        sig_instre.m_overiter(ii,istat) = mean(SITH_TFC(:));
        sig_instre.SD_overiter(ii,istat) = std(SITH_TFC(:));
        sig_instre.iter(ii)=ii*s_step;
        clear SITH_TFC
    end
end
if doplot
figure,
plot(sig_instre.iter,sig_instre.m_overiter,'-','linewidth',2),hold on,
plot(sig_instre.iter,sig_instre.m_overiter+sig_instre.SD_overiter,'--')
plot(sig_instre.iter,sig_instre.m_overiter-sig_instre.SD_overiter,'--')
end
for istat=1:length(statlev)
    SITH_TFC=prctile(cTOT,statlev(istat),3);
    sig_instre.m (istat) = mean(SITH_TFC(:));
    sig_instre.SD (istat) = std(SITH_TFC(:));
    clear SITH_TFC
end
        
        