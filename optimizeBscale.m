%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Gladytz & Joao Periquito
% 20.05.07 (b-value optimization for NNLS & LS simulations)  
%
% This script optimizes the b-values used in the numerical simulations. 
% The b-scale yield a constant signal intensity decrement from one b-value 
% to another. Maintaining the intensity drop constant from one measurement 
% to the next ensures independence of the individual measurements of the 
% signal decay. The functions assumes a tri-exponential decay. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This function optmizes the b-scale in equidistant points  
% - The equidistant intensity drop b-scale is computed by interpolating on 
% a fine evaluation of the tri-exponential decay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables and Descriptions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - decayparamsForOpt - parameters of the tri-exponential. should contain 
% the weights of the exponentials as first column and decay rates as second 
% column
% 
% e.g. [frac_fastOpt diff_fast; frac_medOpt diff_med; frac_slowOpt diff_slow]
% diff_fast   = 0.180;
% diff_med    = 0.0058;
% diff_slow   = 0.0015;
% frac_fastOpt = 0.075;
% frac_medOpt = 0.40;
% frac_slowOpt = 0.525;
%
% - numB - desired number of b-values
%
% - input variables:
%   - decayparamsForOpt
%   - numB
% - output variables:
%   - bscale - optimized b-scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bscale] = optimizeBscale(decayparamsForOpt,numB)
  % [bscale] = optimizeBscale([0.08 0.18; 0.52 0.0058; 0.4 0.0015],30);
  if sum(decayparamsForOpt(:,1))~=1;
    warning('The weights of the exponentials should sum to 1.')
  end
  
  iniScaleEnd = -log(1./numB)./min(decayparamsForOpt(:,2));
  iniScaleStep = -log(1-1./numB)./max(decayparamsForOpt(:,2));
  iniBscale=(0:sqrt(iniScaleStep):sqrt(iniScaleEnd)).^2; % quadratic scale with much to many values for a start.
  intensityScale=1:(-1./numB):(1./numB);
  decay=sum(repmat(decayparamsForOpt(:,1),1,length(iniBscale)).*exp(-decayparamsForOpt(:,2)*iniBscale),1);
  bscale =round(interp1(decay,iniBscale,intensityScale));
  decay2=sum(repmat(decayparamsForOpt(:,1),1,numB).*exp(-decayparamsForOpt(:,2)*bscale),1);
end
