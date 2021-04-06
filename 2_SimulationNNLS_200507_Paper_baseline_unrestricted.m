%% Simulation - Steps: 300; NumSimulations =500; Baseline parameters; unrestricted b-value
load noise_matrix.mat noise
%noise_matrix

%% Steps 300
tic
smallreal=double(realmin('single')).^2;

% numOfSimulations = 100;
number_of_b_valuesDONE=0; %variable used to check the progress of the simulation ----------TG: should be 0 in the beginning
%% Diffusion Parameters (Baseline)
diff_fast   = 0.180;
diff_med    = 0.0058; %[0.006 0.007 0.008 0.009 0.010];
diff_slow   = 0.0015;
frac_fast   = 0.10;
frac_med    = 0.30;
frac_slow   = 0.60;

frac_fastOpt = 0.075;
frac_medOpt = 0.40;
frac_slowOpt = 0.525;

ADCThresh=1./sqrt([diff_fast*diff_med diff_med*diff_slow]);
%% Generate NNLS space of values
ADCBasisSteps = 300;
ADCBasis = logspace( log10(5), log10(2200), ADCBasisSteps);

ADCBasis_fix = [1/diff_fast logspace( log10(ADCThresh(1)), log10(ADCThresh(2)), ADCBasisSteps/3-2) 1/diff_slow];


%% b-value selection
% number_of_b_values =10:5:50;

%% SNR selection 
SNR = [40 80 120 160 200 280 360 440 520 640 760 880 1000];

%% Needed variables creation
ImgTest = zeros(length(number_of_b_values),numOfSimulations,size(SNR,2),max(number_of_b_values)); % TG: One could switch 2nd and 3rd dimension to be consistent
amplitudes = zeros(length(number_of_b_values),numOfSimulations, size(SNR,2), ADCBasisSteps);
amplitudes_fix = zeros(length(number_of_b_values),numOfSimulations, size(SNR,2), ADCBasisSteps/2);
resnorm = zeros(length(number_of_b_values), numOfSimulations, size(SNR,2),1);
resnorm_fix = zeros(length(number_of_b_values), numOfSimulations, size(SNR,2),1);
y_recon = zeros(length(number_of_b_values), numOfSimulations, size(SNR,2), max(number_of_b_values));
y_recon_fix = zeros(length(number_of_b_values), numOfSimulations, size(SNR,2), max(number_of_b_values));
resid = zeros(length(number_of_b_values),numOfSimulations, size(SNR,2), length(number_of_b_values));
resid_fix = zeros(length(number_of_b_values),numOfSimulations, size(SNR,2), length(number_of_b_values));
resultsPeaks = zeros(9,length(number_of_b_values), size(SNR,2), numOfSimulations);
resultsPeaks_fix = zeros(9,length(number_of_b_values), size(SNR,2), numOfSimulations);
HeatMap = zeros (7,length(number_of_b_values),size(SNR,2));
list_of_b_values = zeros(length(number_of_b_values),max(number_of_b_values));

for h=1:length(number_of_b_values)

    
%% b-value optimization  
[b_values] = optimizeBscale([frac_fastOpt diff_fast; frac_medOpt diff_med; frac_slowOpt diff_slow],number_of_b_values(h));

list_of_b_values(h,1:length(b_values)) = b_values;

%% Generate NNLS space of values
A = exp( -kron(b_values',1./ADCBasis));
A_fix = exp( -kron(b_values',1./ADCBasis_fix));


%% Generate Tri-exponential
SI =  ((1-frac_med-frac_fast)*exp(-b_values*diff_slow)+frac_med*exp(-b_values*diff_med)+frac_fast*exp(-b_values*diff_fast));
noiseLevel = SI(1)./SNR;
    for k=1:numOfSimulations
    %% Create Syntetic Data: Perfect data + rician noise
        for i=1:size(noiseLevel,2)
            %ImgTest(h,i,k,1:length(SI)) = SI + noiseLevel(i) * randn(size(SI)); % TG: This is gaussian noise. Shouldn't rician noise have a contribution perpendicular to the true signal? 
            ImgTest(h,k,i,1:length(SI)) = abs(SI + noiseLevel(i) * noise(k,1:number_of_b_values(h)));
            % Rician Bias Correction
            noiseLevEstim=noiseLevel(i)*abs(1+0.05*randn(1)); % 5% uncertainty on noiseEstimation corresponds to determining the noise from 100 noisy points
            ImgTest(h,k,i,1:length(SI)) = sqrt(max(ImgTest(h,k,i,1:length(SI)).^2-noiseLevEstim.^2,smallreal));
        end
    end
toc

    for m=1:numOfSimulations
    
        for j=1:size(noiseLevel,2)

             TempImgTest = squeeze(ImgTest(h,m,j,:));
             TempImgTest = TempImgTest(TempImgTest~=0);
            [ TempAmplitudes, TempResnorm, TempResid ] = CVNNLS(A, TempImgTest);
            amplitudes(h,m,j,:) = TempAmplitudes';
            resnorm(h,m,j,:)    = TempResnorm';
            resid(h,m,j,1:length(TempResid))  = TempResid';
            y_recon (h,m,j,1:size(A,1)) = A * TempAmplitudes;
            
            try %check if we detect three peaks and retrieve the value
            [ GeoMeanRegionADC_1,GeoMeanRegionADC_2,GeoMeanRegionADC_3,RegionFraction1,RegionFraction2,RegionFraction3 ] = NNLS_result_mod( TempAmplitudes, ADCBasis);
            resultsPeaks(1,h,j,m) = RegionFraction1; %(frac_fast - RegionFraction1)./frac_fast.*100;
            resultsPeaks(2,h,j,m) = RegionFraction2; %(frac_med - RegionFraction2)./frac_med.*100;
            resultsPeaks(3,h,j,m) = RegionFraction3; %(frac_slow - )./frac_slow.*100;
            resultsPeaks(4,h,j,m) = GeoMeanRegionADC_1; %(diff_fast - GeoMeanRegionADC_1./1000)./diff_fast.*100;
            resultsPeaks(5,h,j,m) = GeoMeanRegionADC_2; %(diff_med - GeoMeanRegionADC_2./1000)./diff_med.*100;
            resultsPeaks(6,h,j,m) = GeoMeanRegionADC_3; %(diff_slow - GeoMeanRegionADC_3./1000)./diff_slow.*100;
            catch
            % TG: I think we should evaluate also the results, where we have more or less than 3 peaks to be fair. But this could be done after the Simulation is finished.
            % One could do something like: 1. find the two widest minima; 2. use the centers of these minima as borders between the components.
            resultsPeaks(1,h,j,m) = NaN;
            resultsPeaks(2,h,j,m) = NaN;
            resultsPeaks(3,h,j,m) = NaN;
            resultsPeaks(4,h,j,m) = NaN;
            resultsPeaks(5,h,j,m) = NaN;
            resultsPeaks(6,h,j,m) = NaN;
            end
          
            [ TempAmplitudes_fix, TempResnorm_fix, TempResid_fix ] = CVNNLS_PartialRegularization(A_fix, TempImgTest, [1 length(ADCBasis_fix)]);
            amplitudes_fix(h,m,j,:) = TempAmplitudes_fix';
            resnorm_fix(h,m,j,:)    = TempResnorm_fix';
            resid_fix(h,m,j,1:length(TempResid))  = TempResid_fix';
            y_recon_fix (h,m,j,1:size(A,1)) = A_fix * TempAmplitudes_fix;
            
            try %check if we detect three peaks and retrieve the value
            [ GeoMeanRegionADC_1_fix,GeoMeanRegionADC_2_fix,GeoMeanRegionADC_3_fix,RegionFraction1_fix,RegionFraction2_fix,RegionFraction3_fix ] = NNLS_resultTG( TempAmplitudes_fix, ADCBasis_fix ,ADCThresh);
            resultsPeaks_fix(1,h,j,m) = RegionFraction1_fix; %(frac_fast - RegionFraction1)./frac_fast.*100;
            resultsPeaks_fix(2,h,j,m) = RegionFraction2_fix; %(frac_med - RegionFraction2)./frac_med.*100;
            resultsPeaks_fix(3,h,j,m) = RegionFraction3_fix; %(frac_slow - )./frac_slow.*100;
            resultsPeaks_fix(4,h,j,m) = GeoMeanRegionADC_1_fix; %(diff_fast - GeoMeanRegionADC_1./1000)./diff_fast.*100;
            resultsPeaks_fix(5,h,j,m) = GeoMeanRegionADC_2_fix; %(diff_med - GeoMeanRegionADC_2./1000)./diff_med.*100;
            resultsPeaks_fix(6,h,j,m) = GeoMeanRegionADC_3_fix; %(diff_slow - GeoMeanRegionADC_3./1000)./diff_slow.*100;
            catch
            % TG: I think we should evaluate also the results, where we have more or less than 3 peaks to be fair. But this could be done after the Simulation is finished.
            % One could do something like: 1. find the two widest minima; 2. use the centers of these minima as borders between the components.
            resultsPeaks(1,h,j,m) = NaN;
            resultsPeaks(2,h,j,m) = NaN;
            resultsPeaks(3,h,j,m) = NaN;
            resultsPeaks(4,h,j,m) = NaN;
            resultsPeaks(5,h,j,m) = NaN;
            resultsPeaks(6,h,j,m) = NaN;
            end

        end
    end
toc    
    %% display how many selection of b-values are needed
    number_of_b_valuesDONE =number_of_b_valuesDONE+1;
    length(number_of_b_values)-number_of_b_valuesDONE
    save(['steps300_sim500_baseline_unres' datestr(now) '.mat'])
    %toc./60
end


for o=1:length(number_of_b_values)
    for n=1:size(noiseLevel,2)
%     keyboard
    HeatMap(1,o,n) = (sum(isnan(resultsPeaks(1,o,n,:)))./numOfSimulations).*100;
    HeatMap(2,o,n) = mean(squeeze(abs(resultsPeaks(1,o,n,:))),'omitnan');
    HeatMap(3,o,n) = mean(squeeze(abs(resultsPeaks(2,o,n,:))),'omitnan');
    HeatMap(4,o,n) = mean(squeeze(abs(resultsPeaks(3,o,n,:))),'omitnan');
    HeatMap(5,o,n) = mean(squeeze(abs(resultsPeaks(4,o,n,:))),'omitnan');
    HeatMap(6,o,n) = mean(squeeze(abs(resultsPeaks(5,o,n,:))),'omitnan');
    HeatMap(7,o,n) = mean(squeeze(abs(resultsPeaks(6,o,n,:))),'omitnan');

    end
end
% figure
% boxplot([squeeze(resultsPeaks(1,1,:)), squeeze(resultsPeaks(1,2,:)), squeeze(resultsPeaks(1,3,:)),squeeze(resultsPeaks(1,4,:)), squeeze(resultsPeaks(1,5,:)),squeeze(resultsPeaks(1,6,:)),squeeze(resultsPeaks(1,7,:)),squeeze(resultsPeaks(1,8,:)),squeeze(resultsPeaks(1,9,:)),squeeze(resultsPeaks(1,10,:)),squeeze(resultsPeaks(1,11,:))], 'Labels',{'20' '40' '60' '80' '100' '130' '170' '250' '400' '650' '1000'});

timeElapsed = toc./60
% figure, imagesc(AAA)
save('steps300_sim500_baseline_unres.mat')



