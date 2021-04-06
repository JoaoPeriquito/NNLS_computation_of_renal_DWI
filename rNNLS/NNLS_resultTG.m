function [ GeoMeanRegionADC_1,GeoMeanRegionADC_2,GeoMeanRegionADC_3,RegionFraction1,RegionFraction2,RegionFraction3 ] = NNLS_resultTG( TempAmplitudes, ADCBasis , ADCThresh)
% finds and integrates Peaks below and above two defined thresholds
ADCThresh=sort(ADCThresh);
ThreshInd(1)=find(ADCBasis>=ADCThresh(1),1);
ThreshInd(2)=find(ADCBasis>=ADCThresh(2),1);
[locsMax, pksMax]=peakseekTG(TempAmplitudes,1,eps);
[locsMin, pksMin]=peakseek(-TempAmplitudes-(min(-TempAmplitudes)));

peaksMax = locsMax(pksMax ~= 0);
pksMax = pksMax(pksMax ~= 0);
%while length(peaksMax) > 3.5 % fuse closest peaks to reduce number to 3.
%	[neglPeakPos,neglPeak]=min(diff(peaksMax)); % find the two closest peaks
%	if pksMax(neglPeak)>pksMax(neglPeak+1) % which of both is bigger?
%		neglPeak=neglPeak+1; % neglect the latter one
%	end
%	peaksMax=peaksMax(setdiff(1:end,neglPeak));
%	pksMax=pksMax(setdiff(1:end,neglPeak));
%end
%if length(peaksMax) == 3 

TotalArea = sum(TempAmplitudes);
%Peak1
  %[junk,locMinID]=min(abs(locsMin-(peaksMax(1)+peaksMax(2))/2));
	%Peak1End=locsMin(locMinID);
  if min(peaksMax)<ThreshInd(1)
    peakSel=find(peaksMax<=ThreshInd(1)); %find all peaks below Threshhold 1;
    [junk,peakSelID]=max(pksMax(peakSel)); % find the highest of those
    Peak1ID=peakSel(peakSelID(1));
    Peak1End = locsMin(locsMin > peaksMax(Peak1ID));
  else % emergency case: Just use all amplitudes below Threshhold 1
    Peak1End = ThreshInd(1);
  end
	range1 = 1:Peak1End(1);
  RegionFraction1 = sum(TempAmplitudes(range1)) / TotalArea;
  GeoMeanRegionADC_1 = 1000./exp(dot(log(ADCBasis(range1)),TempAmplitudes(range1)) ./ sum(TempAmplitudes(range1)));
  
%		Peak1End = locsMin(locsMin > (peaksMax(1)));
%		Peak2End = locsMin(locsMin > (peaksMax(2)));
%        range1 = ADCBasis >= ADCBasis(1)  &  ADCBasis < ADCBasis(Peak1End(1)+1);
%        ADCBasisRange1 = ADCBasis(range1);
%        ADCampsRange1 = TempAmplitudes(range1);
%        RegionFraction1 = sum( ADCampsRange1 ) / TotalArea;
%        ADCwidth1 = ADCBasisRange1(ADCampsRange1 >= pksMax(1)./2);
%        ADCwidth1 = 1./ADCwidth1(1) - 1./ADCwidth1(end);
%        GeoMeanRegionADC_1 = (1./ exp( dot( ADCampsRange1, log( ADCBasisRange1 ) ) ./ ( RegionFraction1*TotalArea ) )).*1000;

%Peak2
  %[junk,locMinID]=min(abs(locsMin-(peaksMax(2)+peaksMax(3))/2));
	%Peak2End=locsMin(locMinID):
  peakSel=find(peaksMax>=ThreshInd(1)&peaksMax<ThreshInd(2)); %find all peaks between Threshhold 1 and 2;
  if ~isempty(peakSel)
    [junk,peakSelID]=max(pksMax(peakSel)); % find the highest of those
    Peak2ID=peakSel(peakSelID(1));
    Peak2Beg = locsMin(locsMin < peaksMax(Peak2ID));
    Peak2End = locsMin(locsMin > peaksMax(Peak2ID));
  else % emergency case: Just use remaining amplitudes below Threshhold 2
    Peak2Beg = Peak1End+1;
    Peak2End = ThreshInd(2);
  end
	range2 = Peak2Beg(end):Peak2End(1);
  RegionFraction2 = sum(TempAmplitudes(range2)) / TotalArea;
  GeoMeanRegionADC_2 = 1000./exp(dot(log(ADCBasis(range2)),TempAmplitudes(range2)) ./ sum(TempAmplitudes(range2)));
	
%        range2 = ADCBasis >= ADCBasis(Peak1End(1)+1)  &  ADCBasis < ADCBasis(Peak2End(1)+1);
%        ADCBasisRange2 = ADCBasis(range2);
%        ADCampsRange2 = TempAmplitudes(range2);
%        RegionFraction2 = sum( ADCampsRange2 ) / TotalArea;
%        ADCwidth2 = ADCBasisRange2(ADCampsRange2 >= pksMax(2)./2);
%        ADCwidth2 = 1./ADCwidth2(1) - 1./ADCwidth2(end);
%        GeoMeanRegionADC_2 = (1./exp( dot( ADCampsRange2, log( ADCBasisRange2 ) ) ./ ( RegionFraction2*TotalArea ) )).*1000;

 %Peak3
  peakSel=find(peaksMax>=ThreshInd(2)); %find all peaks between Threshhold 1 and 2;
  if ~isempty(peakSel)
    [junk,peakSelID]=max(pksMax(peakSel)); % find the highest of those
    Peak3ID=peakSel(peakSelID(1));
    Peak3Beg = locsMin(locsMin < peaksMax(Peak3ID));
  else % emergency case: Just use all remaining amplitudes
    Peak3Beg = Peak2End+1;
  end
	range3 = Peak3Beg(end):length(ADCBasis);
  RegionFraction3 = sum(TempAmplitudes(range3)) / TotalArea;
  GeoMeanRegionADC_3 = 1000./exp(dot(log(ADCBasis(range3)),TempAmplitudes(range3)) ./ sum(TempAmplitudes(range3)));

%        range3 = ADCBasis >= ADCBasis(Peak2End(1)+1)  &  ADCBasis < ADCBasis(end);
%        ADCBasisRange3 = ADCBasis(range3);
%        ADCampsRange3 = TempAmplitudes(range3);
%        RegionFraction3 = sum( ADCampsRange3 ) / TotalArea;
%        GeoMeanRegionADC_3 = (1./exp( dot( ADCampsRange3, log( ADCBasisRange3 ) ) ./ ( RegionFraction3*TotalArea ) )).*1000;
%else
%    resultsPeaks(1) = 0;
%end

end