%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Gladytz & Joao Periquito
% 20.05.07 (Noise Matrix creation for NNLS & LS simulations)  
%
% This script creates the space of values used for NNLS and LS simulations
% In order to keep the comparison between different SNR-levels and simulated
% scenarios free from statistical fluctuations, the same random numbers were 
% used for all SNR-levels and scenarios. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This function generates the necessary random numbers.
% - Dimensions [number of simulations, maximum number of b-values]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables and Descriptions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - G is used interally and can be used to check how lambda affects G. The
% CV approach attempts to minimize G using lambda.
%
% CV approach based on Golub et al. Technometrics:21; 215-23 (1979).
% Cross-Validation as a Method for Choosing a Good Ridge Parameter
%
% Bisection method used to find minimum
% http://en.wikipedia.org/wiki/Bisection_method
%
% - input variables:
%   - numOfSimulations
%   - number_of_b_values
% - output variables:
%   - noise - noise matrix used for all simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

numOfSimulations = 500;
number_of_b_values =10:5:50;

noise=randn(numOfSimulations,max(number_of_b_values)) + sqrt(-1).*randn(numOfSimulations,max(number_of_b_values));

save(['noise_matrix' datestr(now) '.mat']); % just to make sure that the noise Matrix used in the simulations is not accidentally overwritten
