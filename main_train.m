clc;
close all;
clear all;
warning off;


dataFold = uigetdir();
dataOut = dir([dataFold,'\','*.mat']);
desFreq = [zeros 0.4 0.4 ones];
desirMag = [ones ones zeros zeros];

outData = [];


for i = 1:length(dataOut)
    dataSet = dataOut(i).name;
    inData =load([dataFold '\' dataSet]);
    outvalData = inData.val;
    windowSize = 5;
    noofChannel = [size(outvalData,1)-1]/2;
    fprintf(outvalData);
    fs = 20000;
    f0 = 50;
    fn = fs/2;
    freqRatio = f0/fn;
    thersh = 0.13;
    
    notchWidth = 0.1;
    zeros = [exp( sqrt(-1)*pi*freqRatio ),exp( -sqrt(-1)*pi*freqRatio )];
    poles = (1-notchWidth) * zeros;
    
    % Multiplier Design Scenario
    Fs = f0*2;
    for ii = 1
        filterSignal = filter(ones(1,windowSize)/windowSize,1,outvalData(ii,:));
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
        relatInp = ['numofic '];
        [icaextractSig,meanicaSig,maxicaSig,field_coeff] = independent_analys(ratLev',dataOut,relatInp);
        
        %four channel deployment with Eignen reductionlaity
        for b = 1:length(eig_vec)-(size(eig_vec,1)-windowSize*2)
            for q = 1: length(eig_vec)-(size(eig_vec,1)-windowSize*2)
                channel_data{b,q}  = (abs(eig_vec(b,q)));
            end
        end
        freq_channel = channel_data;
        
        %principle component analysis
        [coeff_matrix,laent_output] = pca(icaextractSig');
        laent_output = abs(laent_output);
        Hpsd = dspdata.psd(laent_output','Fs',Fs);
        power = Hpsd.Data;
        div_space = linspace(0.1,1,length(power));
        collectData = power.*div_space';
        single_comp = mean(collectData);
        
        %                 statistical parameters
        mean_data = abs(mean((icaextractSig)))/1e02
        variance_data = abs(var((icaextractSig)))/1e03
        kurtosis_data = abs(kurtosis((icaextractSig)))
        skewness_data = abs(skewness((icaextractSig)))
        
        %statistical parameters
        mean_data_pca = abs(mean((laent_output)))
        variance_data_pca = abs(var((laent_output)))/1e03
        kurtosis_data_pca = abs(kurtosis((laent_output)))
        skewness_data_pca = abs(skewness((laent_output)))
        
        
        % auto regression method
        [Pxx,F] = pyulear(icaextractSig,2,length(laent_output),1);
        fprintf('Principle componenet data is: %0.5f \n',single_comp)
        plot(F,10*log10(Pxx))
        hold on
        drawnow;
        
        
        %Feature Extractions
        mean_deviation= mean(Pxx)/1e05
        median_deviation = std(Pxx)/1e05
        skewness_eval= skewness(Pxx)
        kurtosis_data = kurtosis(Pxx)
    end
    outData = [outData Pxx];
end
outData = outData/1e05;

%Signal Representation
figure,
plot(outvalData(1,:))
xlabel('Frequecy (hz)')
ylabel('Amplitude')
title('EEG signal')

figure,
plot(smooth(new_soludata))
title('Multiplier Reduction Signal')
xlabel('Frequecy (hz)')
ylabel('Amplitude')

figure,
plot(smooth(full_feature))
title('Discrete Wavelet Signal')
xlabel('Frequecy (hz)')
ylabel('Amplitude')
xlim([3 130])


figure,
plot(smooth(ratLev))
title('Time Domain Signal Analysis')
xlabel('Time (Ms)')
ylabel('Amplitude')


figure,
plot(smooth(icaextractSig))
title('ICA Signal Analysis')
xlabel('Time (Ms)')
ylabel('Amplitude')

figure,
plot(smooth(collectData))
title('Power Spectral Analysis with PCA')
xlabel('Time (Ms)')
ylabel('Amplitude')


% xlabel('Frequency (Hz)')
% ylabel('Power 10*log10({\gamma}V2/Hz')
% title('Multiplier + PCA analysis')

rangeData = sort(collectData,'descend');
figure,plot((fliplr(smooth(collectData))))

% xlim([1 4])
% xlabel('Frequency (Hz)')
% ylabel('Power 10*log10({\gamma}V2/Hz')

figure,plot((fliplr(smooth(collectData))))
xlim([4 8])
xlabel('Frequency (Hz)')
ylabel('Power 10*log10({\gamma}V2/Hz')

figure,plot((fliplr(smooth(collectData))))
xlim([8 13])
xlabel('Frequency (Hz)')
ylabel('Power 10*log10({\gamma}V2/Hz')

figure,plot((fliplr(smooth(collectData))))
xlim([13 30])
xlabel('Frequency (Hz)')
ylabel('Power 10*log10({\gamma}V2/Hz')


%Classification part
for l = 1:size(dataOut)
    if mean(outData(:,l)) > thersh
        meanData(l,:) = mean(outData(:,l));
        result_set(l,:) = true;
    else mean(outData(:,l)) < thersh;
        meanData(l,:) = mean(outData(:,l));
        result_set(l,:) = false;
    end
end
[svmStruct,svIndex] = svmtrain(outData,result_set,'showplot',true);
helpdlg('Training is Completed')

save database.mat svmStruct outData result_set meanData svIndex field_coeff dataOut
