function [SITH,paramOUT,cTOT] = SITH_spwvd(varargin);


paramOUT.T = 816*4;    
paramOUT.Fbins = 2048; 
paramOUT.v0 = 0.045;   
paramOUT.tau0 = 0.05;  
paramOUT.lambda = 0.3; 
paramOUT.Ns = 100;     
paramOUT.fs = 4;       
paramOUT.statlev = [90,95,99];   
paramOUT.s_step = 5;    
paramOUT.doplot = 0;
paramOUT.randseed = [];


pOUT_names = fieldnames(paramOUT);

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
            paramOUT.(pOUT_names{keyN}) = varargin{1}.(char(pIN_names(J)));
        end
        J = J + 1;
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
            paramOUT.(pOUT_names{keyN}) = varargin{J+1};
        end
        J = J + 2;
    end
    for id=1:length(pOUT_names)
        eval([char(pOUT_names(id)),'=paramOUT.',char(pOUT_names(id)),';'])
    end
end

[fwhm_t,fwhm_f] = TFspwv_res([1:T],Fbins,v0,tau0,lambda,90,fs);
display(['using lambda=',num2str(lambda),'; v0=',num2str(v0),' (Dt=',num2str(fwhm_t,3),'sec); tau0=',num2str(tau0),' (Df=',num2str(fwhm_f,3),'Hz)']);
display(['#surr =',num2str(Ns),' Fbins=',num2str(Fbins),' T=',num2str(T)]);
Dt = round(fwhm_t*fs);
Df = round(fwhm_f*Fbins/fs);
cTOT = nan(length(Df:Df:Fbins/2),length(Dt:Dt:T-Dt),Ns);
for is = 1:Ns

if exist('randseed')
    rs = randseed*Ns;
else
rs = randi([1,5000],1);
end
s1=RandStream('mrg32k3a','seed',rs+is);
s2=RandStream('mrg32k3a','seed',rs*Ns+is);

clear c2 c
% White noises
x1=hilbert(randn(s1,T,1));
x2=hilbert(randn(s2,T,1));
% Estimating SPWV-TFC
[c2]=TFCspwv(x1,x2,Fbins,v0,tau0,lambda,fs);
out(is,1)= sum(c2(:)>1);out(is,2)= sum(c2(:)<0);
% Reduce dimension
c = sqrt(c2(Df:Df:Fbins/2,Dt:Dt:end-Dt));

display(['iter=',num2str(is),' mean TFC=',num2str(mean(c(:)))]);
cTOT(:,:,is) =c;
end

for ii=1:round(size(cTOT,3)/s_step)
    for istat=1:length(statlev)
        SITH_TFC=prctile(cTOT(:,:,1:ii*s_step),statlev(istat),3);
        SITH.m_overiter(ii,istat) = mean(SITH_TFC(:));
        SITH.SD_overiter(ii,istat) = std(SITH_TFC(:));
        SITH.iter(ii)=ii*s_step;
        clear SITH_TFC
    end
end
if doplot
figure,
plot(SITH.iter,SITH.m_overiter,'-','linewidth',2),hold on,
plot(SITH.iter,SITH.m_overiter+SITH.SD_overiter,'--')
plot(SITH.iter,SITH.m_overiter-SITH.SD_overiter,'--')
end
for istat=1:length(statlev)
    SITH_TFC=prctile(cTOT,statlev(istat),3);
    SITH.m (istat) = mean(SITH_TFC(:));
    SITH.SD (istat) = std(SITH_TFC(:));
    clear SITH_TFC
end
SITH.out = sum(out>0);        
paramOUT.fwhm_t = fwhm_t;
paramOUT.fwhm_f = fwhm_f;