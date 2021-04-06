function [ GeoMeanRegionADC_1,GeoMeanRegionADC_2,GeoMeanRegionADC_3,RegionFraction1,RegionFraction2,RegionFraction3 ] = NNLS_result_mod( TempAmplitudes, ADCBasis )

[locsMax, pksMax]=peakseek(TempAmplitudes,1,realmin);
[locsMin, pksMin]=peakseek(-TempAmplitudes-(min(-TempAmplitudes)));

peaksMax = locsMax(pksMax ~= 0);
pksMax = pksMax(pksMax ~= 0);
if length(peaksMax) < 2.5 % try to find peaks in curvature
  [locsMax, pksMax]=peakseek(-diff(diff(TempAmplitudes)),1,realmin);
  peaksMax = locsMax(pksMax ~= 0);
  pksMax = TempAmplitudes(peaksMax);
end

while length(peaksMax) > 3.5 % fuse closest peaks to reduce number to 3.
	[neglPeakPos,neglPeak]=min(diff(peaksMax(1:end-1))); % find the two closest peaks
	if pksMax(neglPeak)>pksMax(neglPeak+1) % which of both is bigger?
		neglPeak=neglPeak+1; % neglect the latter one
	end
	peaksMax=peaksMax(setdiff(1:end,neglPeak));
	pksMax=pksMax(setdiff(1:end,neglPeak));
end

%Peak1

Peak1End = locsMin(locsMin > peaksMax(1));
Peak2End = locsMin(locsMin > peaksMax(2));
if isempty(Peak2End)
  [locsMinCurv, pksMinCurv]=peakseek(diff(diff(TempAmplitudes)));% Try to find mimimum in Curvature
  Peak2End = locsMinCurv(locsMinCurv > peaksMax(2));
end

TotalArea = sum(TempAmplitudes);

        range1 = 1:Peak1End(1); %ADCBasis >= ADCBasis(1)  &  ADCBasis < ADCBasis(Peak1End(1)+1);
        
        ADCBasisRange1 = ADCBasis(range1);
        ADCampsRange1 = TempAmplitudes(range1);
        RegionFraction1 = sum( ADCampsRange1 ) / TotalArea;
        
        ADCwidth1 = ADCBasisRange1(ADCampsRange1 >= pksMax(1)./2);
        ADCwidth1 = 1./ADCwidth1(1) - 1./ADCwidth1(end);
        
        GeoMeanRegionADC_1 = (1./ exp( dot( ADCampsRange1, log( ADCBasisRange1 ) ) ./ ( RegionFraction1*TotalArea ) )).*1000;

 %Peak2
        range2 = (Peak1End(1)+1):Peak2End(1); %ADCBasis >= ADCBasis(Peak1End(1)+1)  &  ADCBasis < ADCBasis(Peak2End(1)+1);
        
        ADCBasisRange2 = ADCBasis(range2);
        ADCampsRange2 = TempAmplitudes(range2);
        RegionFraction2 = sum( ADCampsRange2 ) / TotalArea;
        
        ADCwidth2 = ADCBasisRange2(ADCampsRange2 >= pksMax(2)./2);
        if isempty(ADCwidth2)

        else
          ADCwidth2 = 1./ADCwidth2(1) - 1./ADCwidth2(end);
        end;
        
        GeoMeanRegionADC_2 = (1./exp( dot( ADCampsRange2, log( ADCBasisRange2 ) ) ./ ( RegionFraction2*TotalArea ) )).*1000;

 %Peak3
        range3 = (Peak2End(1)+1):length(ADCBasis); %ADCBasis >= ADCBasis(Peak2End(1)+1)  &  ADCBasis < ADCBasis(end);
        
        ADCBasisRange3 = ADCBasis(range3);
        ADCampsRange3 = TempAmplitudes(range3);
        RegionFraction3 = sum( ADCampsRange3 ) / TotalArea;
        
        GeoMeanRegionADC_3 = (1./exp( dot( ADCampsRange3, log( ADCBasisRange3 ) ) ./ ( RegionFraction3*TotalArea ) )).*1000;
end