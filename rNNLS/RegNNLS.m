%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thorarin Bjarnason
% 2006.08.01 (RegNNLS style)  Last modified 2008.04.07, Confirm variables
% and clean up code
%
% RegNNLS.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This function regularizes the matlab nnls function 'lsqnonneg'.
% - Sample Call <copy and paste somewhere else, making sure this function
%       is in the PATH>:
%{
    %Inital values
    te = [10:10:320];
    decay = 200 * exp( -te/20 ) + 800 * exp( -te/80 ) + 2.4;
    noise_factor = 10;
    decay = decay + 5*randn( size( decay ) );
    %Generate Basis Values
    T2min = te(1)*1.5;
    T2max = 3*te(end);
    T2length = 120;
    T2Times = logspace( log10( T2min ), log10( T2max ), T2length );
    T2Basis = exp( -kron( te',1./T2Times ) );
    %NNLS smoothing info
    ChiMin = 1.02;
    ChiMax = 1.025;
    %Append DC offset
    T2Basis = [ T2Basis ones( length( T2Basis(:,1) ),1 ) ];
    %
    %Run NNLS
    [ amplitudes, resnorm, resid] = ...
        RegNNLS(T2Basis, decay', ChiMin, ChiMax);
    %
    %%See result
    %Remove DC offset
    DCOffset = amplitudes(end);
    amplitudes = amplitudes(1:end-1);
    %Calculate Fit
    y_recon = T2Basis*[amplitudes ; DCOffset ];
    figure(1);
    semilogy( te , decay , 'ko' , te , y_recon , 'b-');
%}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - standard Matlab dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables and Descriptions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input Variables:
% - chi_max = input value. The maximum chi^2 value acceptabel for smoothing
% - chi_min = input value. The minimum chi^2 value acceptable for smoothing
% - NNLS_A = input value. The basis functions that will be fit to the data
% - RawData = input value. The actual raw data being fit
%Return Variables:
% - amplitudes = the resulting answer that is passed to the caller
% - resnorm = residual norm of the fit that is passed to the caller
% - resid = the actual residuals of the fit that is passed to the caller
%Internal Variables:
% - chi.
%     smoothEdge = chi^2_smooth of bisector method's outer edge
%     smoothMid = chi^2_smooth of bisector method's middle
% - chisq = chi^2 of the regularized fit, kicked out as a check
% - exitflag = flag used to exit infinite loop
% - highedge = tmp value indicating the high-edge of mu
% - lowedge = tmp value indicating the low-edge of mu
%    edge = mu value at the edge for bisector method
%     max = initial maximum mu value for bisector method
%     mid = middle mu value for bisector method
%     min = initial minimum mu value for bisector method
% - mu_used = mu used such that chi^2_min is between chi^2_max and 
%       chi^2_min
% - NNLS.
%    A = T2 distribution basis functions used for fitting
%    Areg = regularizing term of the basis functions
%    ExitFlag = lsqnonneg exit flag
%    Resid = residuals of the fit
%    Resnorm = chi^2 value of the fit <including the regularization term>
%    ResnormMin = minimum chi^2 value when fitting without regularization
%    y_recon = fitted curve
%    x = amplitudes of lsqnonneg fit
% - T2.
%    total = length of the T2 distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin function RegNNLS

function [ amplitudes, resnorm, resid ] ...
    = RegNNLS(NNLS_A, RawData, chi_min, chi_max)

%Construct the exponential kernel for NNLS
%Create regularizing portion
%
mu.max = ( RawData(1)-RawData(end) )*0.00001;
mu.min = 0;
NNLS.A = NNLS_A;
clear('NNLS_A');

T2.total = length( NNLS.A(1,:));
%Ao = mu*(-2*eye(K+1,K+1) + diag(ones(K,1),1)+ diag(ones(K,1),-1));
%NNLS.Areg = eye( T2.total );  %area
NNLS.Areg = -2*eye(T2.total,T2.total) + diag(ones(T2.total-1,1),1)+ ...
    diag(ones(T2.total-1,1),-1);   %Curvature
%

%Solve minimum case
[ NNLS.x ] = lsqnonneg(  NNLS.A,RawData  ); 
%determine residuals 
resid = RawData-NNLS.A*NNLS.x;
NNLS.ResnormMin = sum(resid.*resid);
%

chi.smoothEdge = NNLS.ResnormMin;
%Change values to be relative to chi2min
chi_min = NNLS.ResnormMin * chi_min;
chi_max = NNLS.ResnormMin * chi_max;


%Using bisector method to zone in on answer




%determin step size <half way between>
mu.mid = mu.max - (mu.max - mu.min)/2;
%

%Determine initial chi_smooth_mid
[ NNLS.x, BexitFlag] = blocknnls( [NNLS.A;(mu.mid)*NNLS.Areg], ...
       [ RawData;zeros(T2.total,1) ] );
if BexitFlag == 0
       [ NNLS.x ] = lsqnonneg( [NNLS.A;(mu.mid)*NNLS.Areg], ...
           [ RawData;zeros(T2.total,1) ] );
end
NNLS.y_recon = NNLS.A * NNLS.x;
resid = RawData-NNLS.y_recon;
NNLS.Resnorm = sum(resid.*resid);
chi.smoothMid = NNLS.Resnorm;


%Set up inital case, keep walking upwards until chi.smoothMid becomes
%greater than chi_min
while chi.smoothMid < chi_min
    mu.mid = mu.mid*2;
    [ NNLS.x, BexitFlag] = blocknnls( [NNLS.A;(mu.mid)*NNLS.Areg], ...
        [ RawData;zeros(T2.total,1) ] );
    if BexitFlag == 0
        [ NNLS.x ] = lsqnonneg( [NNLS.A;(mu.mid)*NNLS.Areg], ...
            [ RawData;zeros(T2.total,1) ] );
   end
    NNLS.y_recon = NNLS.A * NNLS.x;
    resid = RawData-NNLS.y_recon;
    NNLS.Resnorm = sum(resid.*resid);
    chi.smoothMid = NNLS.Resnorm;
    chisq = chi.smoothMid/NNLS.ResnormMin;  %#ok<NASGU> %uncomment to see
end

exitflag = 0;    
if (chi.smoothMid < chi_max) && (chi.smoothMid > chi_min) %Keeper!
    exitflag = 1;
    mu_used = mu.mid;  %#ok<NASGU> %uncomment to see
    amplitudes = NNLS.x(1:end);
    chisq = chi.smoothMid/NNLS.ResnormMin;
    resnorm = chisq;
    NNLS.y_recon = NNLS.A * NNLS.x;
    resid = NNLS.y_recon - RawData;
else
    mu.edge = mu.mid;
    chi.smoothEdge = chi.smoothMid;
    mu.mid = (mu.max - mu.mid)/2;
end
%

%Loop until convergence
while exitflag == 0
    %Do nnls
    [ NNLS.x, BexitFlag] = blocknnls( [NNLS.A;(mu.mid)*NNLS.Areg], ...
        [ RawData;zeros(T2.total,1) ] );
    if BexitFlag == 0
        [ NNLS.x ] = lsqnonneg( [NNLS.A;(mu.mid)*NNLS.Areg], ...
            [ RawData;zeros(T2.total,1) ] );
   end
    NNLS.y_recon = NNLS.A * NNLS.x;
    resid = RawData-NNLS.y_recon;
    NNLS.Resnorm = sum(resid.*resid);
    chi.smoothMid = NNLS.Resnorm;
    chisq = chi.smoothMid/NNLS.ResnormMin;  %uncomment to see
    %Set cases
    if ( chi.smoothMid >= chi_min ) && ( chi.smoothMid <= chi_max ) %Keeper!
        mu_used = mu.mid;  %#ok<NASGU> %uncomment to see
        amplitudes = NNLS.x(1:end);
        resnorm = chisq;
        resid = NNLS.y_recon - RawData;
        exitflag = 1; %#ok<NASGU>
        break;
    end
    if ( mu.edge > mu.mid ) %edge is upper edge
        if chi.smoothMid > chi_max
            mu.min = mu.edge - 2*(mu.mid-mu.min);
            mu.edge = mu.mid;
            mu.mid = mu.edge - (mu.edge - mu.min)/2;
        else  %chi.smoothMid < chi_min
            mu.max = mu.edge;
            mu.edge = mu.mid;
            mu.mid = mu.max - (mu.max - mu.mid)/2;
        end
    else  %the case mu.edge < mu.mid, so edge is lower edge
        if chi.smoothMid < chi_min
            mu.max = mu.edge + 2*mu.mid;
            mu.edge = mu.mid;
            mu.mid = mu.max - (mu.max - mu.mid)/2;
        else %chi.smooMid > chi_max
            mu.min = mu.edge;
            mu.edge = mu.mid;
            mu.mid = mu.min - (mu.min - mu.mid)/2;
        end
    end
end
%

