%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thorarin Bjarnason
% 2009.01.21 (Regularized NNLS using cross validation)  
% Last modified 2010.05.15
%
% CVNNLS.m
%
% This fitting routine uses fastnnls from the N-way toolbox
% C. A. Andersson and R. Bro. The N-way Toolbox for MATLAB.
%  Chemom.Intell.Lab.Syst. 52 (1):1-4, 2000.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This function regularizes the matlab nnls function 'lsqnonneg'.
% - Sample Call <copy and paste somewhere else, making sure this function
%       is in the PATH>:
%{
    %Inital values
    te = 10:10:320;
    decay = 200 * exp( -te/20 ) + 800 * exp( -te/80 ) + 2.4;
    noise_factor = 10;
    decay = decay + noise_factor*randn( size( decay ) );
    %Generate Basis Values
    T2min = te(1)*1.5;
    T2max = 3*te(end);
    T2length = 120;
    T2Times = logspace( log10( T2min ), log10( T2max ), T2length );
    T2Basis = exp( -kron( te',1./T2Times ) );
    %Run NNLS
    [ amplitudes, mu ] = CVNNLS(T2Basis, decay');
    %
    %%See result
    %Calculate Fit
    y_recon = T2Basis*amplitudes;
    figure(1)
    semilogy( te , decay , 'ko' , te , y_recon , 'b-')
    figure(2)
    semilogx(T2Times,amplitudes)
%}
%
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
%   - A = T2 basis
%   - Decay = multiexponential decay
% - output variables:
%   - xOut = amplitudes of the T2 distribution
%   - ChiUsed = chi^2_smooth/chi^2_min - the ratio of the chi^2 to the
%   minimum chi^2 based on the smoothing parameter solved for using this
%   routine.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin function CVNNLS
function [ xOut, ChiUsed, resid ] = CVNNLS(A,Decay)


%I should be able to define H in here as well.
%Cite cross validation method in header info
% G can be used for plotting purposes. TT = sortrows(G);
% Assuming A is of the form <length_decay,length_T2dist>

%Identity matrix
I = eye(length(Decay));
%Curvature
T2length = length( A(1,:) );
H = -2*eye(T2length,T2length) + diag(ones(T2length-1,1),1)+ ...
    diag(ones(T2length-1,1),-1);   %Curvature
%

%Assuming derivatives at these locations have opposite signs. IE the
%positions chosen bracket the solution
% The ranges specified here (0.00001->8, with tolerance 0.0001) behave well
% for SNR 50 -> 1000. <in simulations defined as the intensity of the first
% echo divided by the standard deviation of the noise>. These ranges could
% be narrowed for more specific SNR. Run some simulations and uncomment the
% debugging code below. 2009.01.21, Thorarin Bjarnason.
LambdaLeftInitial = 0.00001;
LambdaRightInitial = 8;
tol = 0.0001;

LambdaLeft = LambdaLeftInitial;
LambdaRight = LambdaRightInitial;

%LambdaLeft = Decay(end)/Decay(1)/100;
%LambdaRight = Decay(end)/Decay(1)*100
%tol = Decay(end)/Decay(1)/20;

% LambdaLeft = Decay(end)/Decay(1)/10
% LambdaRight = LambdaLeft*3000
% tol = LambdaLeft/2;

%function value at left point
G_Left = getG( A, H, I, LambdaLeft, Decay);
%function at left point + small difference
G_LeftD = getG( A, H, I, LambdaLeft+tol, Decay);
%Determine derivative of Left point
f_left = (G_LeftD - G_Left)/tol;
%Assign to G for plotting
G(1,:) = [ LambdaLeft, G_Left ];


i = 1;
%Find minimum
while( abs( LambdaRight - LambdaLeft) ) > tol
    %determine midpoint of domain
    midpoint = (LambdaRight + LambdaLeft)/2;
    %
    
    %function value at middle point
    %now determine function value
    G_Middle = getG( A, H, I, midpoint, Decay );
    %function at left point + small difference
    G_MiddleD = getG( A, H, I, midpoint+tol, Decay);
    %Determine derivative of Left point
    f_middle = (G_MiddleD - G_Middle)/tol;
    %
    
    %Save G for plotting
    G(i+1,:) = [midpoint,G_Middle];
    %
    
    %A check for debugging. I assume that the original lambda's bracket the
    %minimum. So on the first run, if both derivatives have the same sign,
    %we are in trouble
    if (i == 100) 
        error( [ 'Original choice of Lambda might not bracket minimum. '...
            , 'Debug the code and check G' ] )
        break %#ok<UNRCH>
    end
    %
    
    %Continue with logic
    if f_left*f_middle > 0
        %Throw away left half
        LambdaLeft = midpoint;
        f_left = f_middle;
    else
        %Throw away right half
        LambdaRight = midpoint;
    end
    i = i+1; %increment counter
end
% find solution <note: there might be a way to remove this final
% calculation>
Lambda = midpoint;
xOut = nnlsfit( A, H, Lambda, Decay );
%
% Determin chi2_min and chi2_smooth
% chi2_min
[ amp_min, resnormMin] = lsqnonneg(  A, Decay  ); 
% chi2_smooth
y_recon = A * xOut;
resid = Decay - y_recon;
resnormSmooth = sum(resid.*resid);
ChiUsed = resnormSmooth/resnormMin;




% % Plot G for debugging
% % There are some local mins which become obvious when zoomed in with
% % TestLambda. These local mins genearlly appear within 0.01 of the true
% % minimum
% % First do zoomed in case
% TestLambda = linspace(min(G(:,1)),10*Lambda,75);
% for i = 1:length(TestLambda)
%     TestG(i) = getG(A, H, I, TestLambda(i), Decay);
% end
% G = sortrows(G);
% figure(1111)
% plot(G(:,1),G(:,2), 'b', 'LineWidth', 2)
% hold on
% plot(TestLambda,TestG, 'r--', 'LineWidth', 2)
% plot( [ Lambda Lambda], [min(TestG) max(TestG) ], 'k' )
% hold off
% legend('fit','G','min used')
% % now do zoomed out case
% TestLambda2 = linspace(min(G(:,1)),max(G(:,1)),75);
% for i = 1:length(TestLambda)
%     TestG2(i) = getG(A, H, I, TestLambda2(i), Decay);
% end
% figure(1112)
% plot(G(:,1),G(:,2), 'b', 'LineWidth', 2)
% hold on
% plot(TestLambda2,TestG2, 'r--', 'LineWidth', 2)
% hold off
% legend('fit','G')
% %write values to disk
% Gcurve.Gfit = G;
% Gcurve.GcurveSmall = [TestLambda;TestG];
% Gcurve.GcurveBig = [TestLambda2;TestG2];
% Gcurve.Fit = [[ Lambda Lambda]; [min(TestG) max(TestG)] ];
% save( 'Gcurve.mat', 'Gcurve')
% aaa=1;
%


%Determining lambda function, G
function G = getG( A, H, I, Lambda, Decay )
%Fit decay data
NNLSfit = nnlsfit( A, H, Lambda, Decay );
%Calculate G using CrossValidation method.
G = norm( Decay - A*NNLSfit )^2/...
    trace( I - A*( A'*A + (Lambda)*H'*H )^-1*A' )^2;
%

%NNLS fitting routine
function xOut = nnlsfit( A,H,Lambda,Decay )

[ xOut ] = fastnnls( [A;(Lambda)*H]'*[A;(Lambda)*H], ...
    [A;(Lambda)*H]'*[ Decay;zeros(length(H(:,1)),1) ]);
% [ xOut, BexitFlag] = blocknnls( [A;(Lambda)*H], ...
%        [ Decay;zeros(length(H(:,1)),1) ] );
% if BexitFlag == 0
%        [ xOut ] = lsqnonneg( [A;(Lambda)*H], ...
%            [ Decay;zeros(length(H(:,1)),1) ] );
% end
%
