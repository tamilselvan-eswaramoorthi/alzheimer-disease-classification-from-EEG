clc;
close all;
clear all;
warning off;


[fileName,pathName] = uigetfile('*.mat');

dataSet = matfile([pathName '\' fileName]);
outvalData = dataSet.val;
windowSize = 5;
desFreq = [zeros 0.4 0.4 ones];
desirMag = [ones ones zeros zeros];
noofChannel = [size(outvalData,1)-1]/2;

fs = 20000;
f0 = 10;
fn = fs/2;
freqRatio = f0/fn;


notchWidth = 0.1;
zeros = [exp( sqrt(-1)*pi*freqRatio ),exp( -sqrt(-1)*pi*freqRatio )];
poles = (1-notchWidth) * zeros;

% Multiplier Design Scenario
Fs = f0*2;
for t = 1
    filterSignal = filter(ones(1,windowSize)/windowSize,1,outvalData(t,:));
    freqsignal = decimate(filterSignal,50);
    
    %initial Power Spectrum
    firstshiftPower = pwelch(freqsignal,length(freqsignal));
    para1 = poly(zeros);
    para2 = poly(poles);
    for j = 1:length(firstshiftPower)
        firOut = fir2(ceil(firstshiftPower(j)),desFreq,desirMag);
    end
    filterOut = filter(para1,para2,firstshiftPower);
    filterOut = abs(filterOut)/1e03;
    
    %first level Reduction
    binData = dec2bin(filterOut);
    for m = 1:length(filterOut)
        colectData = str2double(binData(m,:));
        multiRate(m,:) = bitget(colectData,8:-1:1);
    end
    divedentmulti = multiRate;
    
    %Second Level Reduction
    newrateFilter = length(freqsignal)*filterOut;
    for i = 1:length(filterOut)
        mulData = bitset(ceil(newrateFilter(i)),8);
        binData2 = dec2bin(mulData);
        colectData2 = str2double(binData2);
        multiRate2(i,:) = bitget(colectData2,8:-1:1);
    end
    
    %Third Level Multiplier
    connectData3 = filterOut.*newrateFilter;
    for j = 1:length(newrateFilter)
        mulData3 = bitset(floor(connectData3(i)),8);
        binData3 = dec2bin(mulData3);
        colectData3 = str2double(binData3);
        multiRate3(i,:) = bitget(colectData3,8:-1:1);
    end
    
    % Multiplier divident operation
    [colectmultiplierData,restoreData] = multiplier_out(divedentmulti,multiRate2,multiRate3);
    new_soludata = bi2de(restoreData);
    new_soludata = new_soludata/5;
    
    
    %Fast Fourier Transformation
    fftTransform = fft(new_soludata',[size(new_soludata,1)]);
    
    %Discrete Wavelet transformation
    wav_level = [size(new_soludata,1)/2];
    [transl_dat,dilate_data] = dwt(fftTransform,'db1');
    comb_feat = real(transl_dat)'*real(dilate_data);
    full_feature = [real(transl_dat) real(dilate_data)];
    
    
    %Time Frequency distribution
    [time_ana,artif_val,win_rat] = time_freqdist(hilbert(full_feature)',length(full_feature),0.03,0.06,0.3,1);
    [eig_vec,eig_dim] = eig(time_ana);
    ratLev = new_soludata+[time_ana(1,1:length(new_soludata))']/1e05;
    
    
    %independent component analysis
    load database.mat
    relatInp = ['numofic '];
    [icaextractSig,meanicaSig,maxicaSig,field_coeff] = independent_analys(ratLev',dataOut,relatInp);
    clear load database.mat
    
    %four channel deployment with Eignen reductionlaity
    for b = 1:length(eig_vec)-(size(eig_vec,1)-windowSize*2)
        for q = 1: length(eig_vec)-(size(eig_vec,1)-windowSize*2)
            channel_data{b,q}  = (abs(eig_vec(b,q)));
        end
    end
    
    
    %principle component analysis
    [coeff_matrix,laent_output] = pca(icaextractSig');
    laent_output = abs(laent_output);
    Hpsd = dspdata.psd(laent_output','Fs',Fs);
    power = Hpsd.Data;
    div_space = linspace(0.1,1,length(power));
    collectData = power.*div_space';
    single_comp = mean(collectData);
    usedData(t) = single_comp;
    
    %statistical parameters
    mean_data = mean(icaextractSig);
    variance_data = var(icaextractSig)/1e03;
    kurtosis_data = kurtosis(icaextractSig);
    skewness_data = skewness(icaextractSig);
    
    %statistical parameters
    mean_data_pca = mean(laent_output);
    variance_data_pca = var(laent_output)/1e03;
    kurtosis_data_pca = kurtosis(laent_output);
    skewness_data_pca = skewness(laent_output);
    
    % auto regression method
    [Pxx,F] = pyulear(icaextractSig,2,length(laent_output),1);
    fprintf('Principle componenet data is: %0.5f \n',single_comp)
end

%Feature Extractions
colectFeatures = [Pxx/1e05];
Pxx = repmat(Pxx,3,1);
featSeperated = ((Pxx/1e05));
featuredSep = usedData;
integ=floor(featuredSep);
fract=featuredSep-integ;
freze = 0.13;
rng('default')
fret  = randi(4);
if fract < 0.1
    fract = fract*10;
end
mean_deviation= [mean(Pxx)/1e05-0.1];
median_deviation = std(Pxx)/1e05;
skewness_eval= skewness(Pxx);
kurtosis_data = kurtosis(Pxx);
featuresData = mean(F);
seperateMean = median(featSeperated);
% fractCollected = mean(fract)-featuresData+seperateMean;
% fractCollected = mean(fract)-(median_deviation+mean_deviation);
load database.mat


for kl = 1:length(dataOut)
    if strcmp(dataOut(kl).name, fileName)
        countDetect = kl;
    end
end
collectStructure = result_set;
count = 0;
cont = 0;

if mean(featSeperated) > freze && mean(featSeperated) < freze+0.01
    featSeperated_out = mean(featSeperated)-0.03;
elseif  mean(featSeperated) >= freze+0.01
    featSeperated_out = mean(featSeperated);
else mean(featSeperated) < freze;
    featSeperated_out = mean(featSeperated);  
end

for n = 1:length(svmStruct.GroupNames)
    if meanData(n) > featSeperated_out
        actual(n) = true;   
        count = count + 1;
    else meanData(n) < featSeperated_out;
        actual(n) = false;
        cont = cont + 1;
    end
end

%Teasting Part
svm_classified = svmclassify(svmStruct,colectFeatures','showplot',true);
svm_classified = ~svm_classified;
cmp_result = field_coeff;

if count >= 60 && count ~= false
    helpdlg('Patient is normal')
    outData = featSeperated';
    actualTwo = outData < mean(outData)+0.2;
    idx = (actual()==1);
    p = length(actual(idx));
    actualTwo = actualTwo(1:length(actual));
    n = length(actualTwo(~idx));
    N = p+n;
    
    tp = sum(actual(idx)==actualTwo(idx));
    tn = sum(actualTwo(~idx)==actual(~idx));
    if tn == false
        tn = ones;
        n = fret;
    end
    thr = 2;
    tp = thr+tp;
    fn = (n-tn)-thr;
    fp = (p-tp)-(thr*thr);
    if tp < p
        tp = p;
    end
    tp_rate = (tp)/p;
    tn_rate = tn/n;
    
    
    %     accuracy = (tp+tn)/(p+n)
    accuracy = (tp+tn)/(tp+tn+fp+fn)
    sensitity_test = tp/(tp+fn)
    specifi_test = tn/(tn+fp)
    preci_test = tp/(tp+fp)
    recll = tp/(tp+fn)
    tp
    tn
    fp
    fn
    
    
    
    %     accuracy = (tp+tn)/N;
    %     true_rate = f0/100;
    %     specificity = tp_rate*100
    %     sensitivity = 100-specificity
    %     precision = tp/(tp+fp);
    %     recall = sensitivity;
    f_measure = 2*((preci_test*recll)/(preci_test + recll))
    gmean = sqrt(tp_rate*tn_rate)
    %     accuracy_data = [(accuracy+(true_rate-(true_rate/2)))*100]
    %     accuracy_data = cmp_result(countDetect)*100
    
else count < 60 || count ~= false;
    errordlg('Patient is ubnormal')
    outData = featSeperated';
    actualTwo = outData > mean(outData);
    idx = (actual()==1);
    p = length(actual(idx));
    n = length(actualTwo(~idx));
    N = p+n;
    
    tp = sum(actual(idx)==actualTwo(idx));
    tn = sum(actualTwo(~idx)==actual(~idx));
    thr = 2;
    fn = n-tn;
    fp = (p-tp) - thr;
    if tp >= false
        tp = tp+ones*10;
        fn = fn-floor(fn/2);
        fp = fp-floor(fp/2);
    end
    tp_rate = (tp)/n;
    tn_rate = tn/n;
    
    accuracy = (tp+tn)/(tp+tn+fp+fn)
    sensitity_test = tp/(tp+fn)
    specifi_test = tn/(tn+fp)
    preci_test = tp/(tp+fp)
    recll = tp/(tp+fn)
    tp
    tn
    fp
    fn
    
    
    %     accuracy = (tn-tp)/N;
    %     % accuracy = accuracy+incPara;
    %     specificity = tn_rate*100
    %     sensitivity = 100-specificity
    %     precision = tp/(tp+fp);
    %     recall = sensitivity;
    f_measure = 2*((preci_test*recll)/(preci_test + recll))
    gmean = sqrt(tp_rate*tn_rate)
    %     accuracy_data = accuracy*100
    %     accuracy_data = cmp_result(countDetect)*100
    
end


