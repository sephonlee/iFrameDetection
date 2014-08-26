function iFramesInfo = detectChangePointBeta(X, tail, isPlot)


if nargin < 3
    isPlot = 0;
end

addpath(genpath('//psf/Home/Desktop/Research/PNNL/sourceCode'))
% Demonstration of online detection of a change in 1d Gaussian parameters.
%
% Implementation of:
% @TECHREPORT{ adams-mackay-2007,
%    AUTHOR = {Ryan Prescott Adams and David J.C. MacKay},
%    TITLE  = "{B}ayesian Online Changepoint Detection",
%    INSTITUTION = "University of Cambridge",
%    ADDRESS = "Cambridge, UK",
%    YEAR = "2007",
%    NOTE = "arXiv:0710.3742v1 [stat.ML]"
% }
%
% Thanks to Ryan Turner and Miguel Lazaro Gredilla for pointing out bugs
% in this.

% First, we wil specify the prior.  We will then generate some fake data
% from the prior specification.  We will then perform inference. Then
% we'll plot some things.

% Start with a clean slate.
% clear;

% How many time steps to generate?
% T = length;
T = size(X, 1);

% Specify the hazard function.
% This is a handle to a function that takes one argument - the number of
% time increments since the last changepoint - and returns a value in
% the interval [0,1] that is the probability of changepoint.  Generally
% you might want to have your hazard function take parameters, so using
% an anonymous function is helpful.  We're going to just use the simple
% constant-rate hazard function that gives geomtrically-drawn intervals
% between changepoints.  We'll specify the rate via a mean.
lambda        = 200;
hazard_func  = @(r) constant_hazard(r, lambda);

% This data is Gaussian with unknown mean and variance.  We are going to
% use the standard conjugate prior of a normal-inverse-gamma.  Note that
% one cannot use non-informative priors for changepoint detection in
% this construction.  The NIG yields a closed-form predictive
% distribution, which makes it easy to use in this context.  There are
% lots of references out there for doing this kind of inference - for
% example Chris Bishop's "Pattern Recognition and Machine Learning" in
% Chapter 2.  Also, Kevin Murphy's lecture notes.
mu0    = 0;
kappa0 = 1;
alpha0 = 1;
beta0  = 1;


% Store the times of changepoints.  It's useful to see them.
CP = [0];


% Now we have some data in X and it's time to perform inference.
% First, setup the matrix that will hold our beliefs about the current
% run lengths.  We'll initialize it all to zero at first.  Obviously
% we're assuming here that we know how long we're going to do the
% inference.  You can imagine other data structures that don't make that
% assumption (e.g. linked lists).  We're doing this because it's easy.
R = [];

% At time t=1, we actually have complete knowledge about the run
% length.  It is definitely zero.  See the paper for other possible
% boundary conditions.
R(1,1) = 1;

% Track the current set of parameters.  These start out at the prior and
% accumulate data as we proceed.
muT    = mu0;
kappaT = kappa0;
alphaT = alpha0;
betaT  = beta0;

% Keep track of the maximums.
maxes  = zeros(T+1,1);
importantPoint = [];

P = 1;
currentColumn = 1;
headColumnHeight = 1;
removeColumns = 0;
% tail = 100;
tailLength = 1;
maxSize = 1;
maxColumnSize = 1;
% Loop over the data like we're seeing it all for the first time.
for t=1:T
    
    %     if t~=1
    %         R(:,t-1) = 0;
    %     end
       
    
    % Evaluate the predictive distribution for the new datum under each of
    % the parameters.  This is the standard thing from Bayesian inference.
    predprobs = studentpdf(X(t), muT, ...
        betaT.*(kappaT+1)./(alphaT.*kappaT), ...
        2 * alphaT);
    
    % Evaluate the hazard function for this interval.
    H = hazard_func([1:size(R, 1)]');
    
    % Evaluate the growth probabilities - shift the probabilities down and to
    % the right, scaled by the hazard function and the predictive
    % probabilities.
%     K = R;
%     R = [];
%     R = [R K];
    R(2:size(R,1)+1,currentColumn+1) = R(:,currentColumn) .* predprobs .* (1-H);
    
    % Evaluate the probability that there *was* a changepoint and we're
    % accumulating the mass back down at r = 0.
    R(1,currentColumn+1) = sum( R(1:end-1,currentColumn) .* predprobs .* H );
    
    % Renormalize the run length probabilities for improved numerical
    % stability.
   
    
    if tailLength == tail
        
%         break;
        
        R(2:headColumnHeight, :) = [];
        removeColumns = removeColumns + (headColumnHeight - 2 + 1);
        muT(1:headColumnHeight-1, :) = [];
        kappaT(1:headColumnHeight-1, :) = [];
        alphaT(1:headColumnHeight-1, :) = [];
        betaT(1:headColumnHeight-1, :) = [];
    end
    
    R(:,currentColumn+1) = R(:,currentColumn+1) ./ sum(R(:,currentColumn+1));
    
    
    P = [P; R(end,currentColumn+1)];
    % Update the parameter sets for each possible run length.
    muT0    = [ mu0    ; (kappaT.*muT + X(t)) ./ (kappaT+1) ];
    kappaT0 = [ kappa0 ; kappaT + 1 ];
    alphaT0 = [ alpha0 ; alphaT + 0.5 ];
    betaT0  = [ beta0  ; betaT + (kappaT .*(X(t)-muT).^2)./(2*(kappaT+1)) ];
    muT     = muT0;
    kappaT  = kappaT0;
    alphaT  = alphaT0;
    betaT   = betaT0;
    
    
%     if count == tail   
%         muT(2:2+cutRows,:) = [];
%         kappaT(2:2+cutRows,:) = [];
%         alphaT(2:2+cutRows,:) = [];
%         betaT(2:2+cutRows,:) = [];
%     end
    
    
    % Store the maximum, to plot later.
    maxes(t) = find(R(:,currentColumn)==max(R(:,currentColumn))) + removeColumns;
    
    if t > 1 && maxes(t,1) - maxes(t-1,1) < 1
        % find a change point
        importantPoint = [importantPoint; t];
        
        % cut tail
        R = R(:,currentColumn:end);
        
        currentColumn = 2;
        tailLength = 1;
        headColumnHeight = size(R, 1) - 1;
    else
        currentColumn = currentColumn + 1;
        tailLength = tailLength + 1;
    end
    
    dim = size(R);
    maxSize = max(maxSize, dim(1)*dim(2));
    maxColumnSize = max(maxColumnSize, dim(1));
end

if isPlot
    % Plot the data and we'll have a look.
    figure(1)
    subplot(2,1,1);
    plot([1:T]', X, 'b-', CP, zeros(size(CP)), 'rx');
    grid;
    
    subplot(2,1,1)
    hold on;
    plot (importantPoint, X(importantPoint), '^r');
    
    % Show the log smears and the maximums.
    subplot(2,1,2);
    colormap(gray());
    imagesc(-log(R));
    hold on;
    plot([1:T+1], [maxes zeros(T+1,T)], 'r-');
    hold off;
end

iFramesInfo.P = P;
iFramesInfo.R = R;
iFramesInfo.maxes = maxes;
iFramesInfo.importantPoint = importantPoint;
iFramesInfo.maxSize = maxSize;
iFramesInfo.maxColumnSize = maxColumnSize;

end