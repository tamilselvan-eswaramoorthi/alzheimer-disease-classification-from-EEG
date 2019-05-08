function [degra_signal,parametric_signals] = SDTH_spwvd(var_signal,var_signal2,varargin)

parametric_signals.T = length(var_signal(:));     
parametric_signals.Fbins = 2048;          
parametric_signals.v0 = 0.045;       
parametric_signals.tau0 = 0.05;        
parametric_signals.lambda = 0.3;         
parametric_signals.Ns = 100;             
parametric_signals.fs = 4;                
parametric_signals.statlev = [90,95,99];  
parametric_signals.s_step = 5;           
parametric_signals.doplot = 0;
parametric_signals.t = [1:parametric_signals.T]/parametric_signals.fs; % Time [s]
parametric_signals.F = linspace(0,parametric_signals.fs/2,parametric_signals.Fbins); % Frequency [Hz]

pOUT_names = fieldnames(parametric_signals);

if isempty(varargin)
for id=1:length(pOUT_names)
     eval([char(pOUT_names(id)),'=paramOUT.',char(pOUT_names(id)),';'])
end   
elseif isstruct(varargin{1})
pIN_names = fieldnames(varargin{1});
J=1;
while J <= length(pIN_names)
    keyN = find(strcmpi(pIN_names(J), pOUT_names));
    if isempty(keyN)
        warning('PPLIKEL:UnknownOption', 'Do not know option ''%s''', key);
    else
        parametric_signals.(pOUT_names{keyN}) = varargin{1}.(char(pIN_names(J)));
    end
    J = J + 2;
end
for id=1:length(pOUT_names)
     eval([char(pOUT_names(id)),'=paramOUT.',char(pOUT_names(id)),';'])
end    
    
else
J = 1;
while J <= length(varargin)
    key = varargin{J};
    keyN = find(strcmpi(key, pOUT_names));
    if isempty(keyN)
        warning('PPLIKEL:UnknownOption', 'Do not know option ''%s''', key);
    else
        parametric_signals.(pOUT_names{keyN}) = varargin{J+1};
    end
    J = J + 2;
end
for id=1:length(pOUT_names)
     eval([char(pOUT_names(id)),'=paramOUT.',char(pOUT_names(id)),';'])
end
end

% Estimating time-frequency resolution
[fwhm_t,fwhm_f] = TFspwv_res([1:T],Fbins,v0,tau0,lambda,90,fs);
display(['using lambda=',num2str(lambda),'; v0=',num2str(v0),' (Dt=',num2str(fwhm_t,3),'sec); tau0=',num2str(tau0),' (Df=',num2str(fwhm_f,3),'Hz)']);
display(['#surr =',num2str(Ns),' Fbins=',num2str(Fbins),' T=',num2str(T)]);
Dt = round(fwhm_t*fs);
Df = round(fwhm_f*Fbins/fs);
c_x1e1_TOT = nan(Fbins/(fs/2),T,Ns); % up to 1 Hz
c_x2e2_TOT = nan(Fbins/(fs/2),T,Ns); % up to 1 Hz
for is = 1:Ns

s1=RandStream('mrg32k3a','seed',is);
s2=RandStream('mrg32k3a','seed',is+Ns);
clear c2 c
% White noises
e1=hilbert(randn(s1,T,1));
e2=hilbert(randn(s2,T,1));
% Estimating SPWV-TFC
[c2_x1e1]=TFCspwv(var_signal,e1,Fbins,v0,tau0,lambda,fs);
[c2_x2e2]=TFCspwv(var_signal2,e2,Fbins,v0,tau0,lambda,fs);
% Reduce dimension
c_x1e1_TOT(:,:,is) = sqrt(c2_x1e1(1:Fbins/(fs/2),:));
c_x2e2_TOT(:,:,is) = sqrt(c2_x2e2(1:Fbins/(fs/2),:));
outx1(is,1)= sum(c2_x1e1(:)>1);outx1(is,2)= sum(c2_x1e1(:)<0);
outx2(is,1)= sum(c2_x2e2(:)>1);outx2(is,2)= sum(c2_x2e2(:)<0);
clear c2_x1e1 c2_x2e2

display(['iter=',num2str(is)]);
end

for ii=1:round(size(c_x1e1_TOT,3)/s_step)
    for istat=1:length(statlev)
        SDTH_TFCx1=prctile(c_x1e1_TOT(:,:,1:ii*s_step),statlev(istat),3);
        SDTH_TFCx2=prctile(c_x2e2_TOT(:,:,1:ii*s_step),statlev(istat),3);
        SDTH_TFC = max(SDTH_TFCx1,SDTH_TFCx2);
        degra_signal.m_overiter(ii,istat) = mean(SDTH_TFC(:));
        degra_signal.SD_overiter(ii,istat) = std(SDTH_TFC(:));
        degra_signal.iter(ii)=ii*s_step;
        clear SDTH_TFC
    end
end
if doplot
figure,
plot(degra_signal.iter,degra_signal.m_overiter,'-','linewidth',2),hold on,
plot(degra_signal.iter,degra_signal.m_overiter+degra_signal.SD_overiter,'--')
plot(degra_signal.iter,degra_signal.m_overiter-degra_signal.SD_overiter,'--')
end
for istat=1:length(statlev)
    SDTH_TFCx1=prctile(c_x1e1_TOT,statlev(istat),3);
    SDTH_TFCx2=prctile(c_x2e2_TOT,statlev(istat),3);
    SDTH_TFC = max(SDTH_TFCx1,SDTH_TFCx2);
    degra_signal.m (istat) = mean(SDTH_TFC(:));
    degra_signal.SD (istat) = std(SDTH_TFC(:));
    degra_signal.TFR(:,:,istat) =SDTH_TFC;
    clear SDTH_TFC
end
degra_signal.outx2=outx2;
degra_signal.outx1=outx1;
        
        