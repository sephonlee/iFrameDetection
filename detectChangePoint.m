function iFramesInfo = detectChangePoint(X, isPlot)


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


% Plot the data and we'll have a look.
figure(2)
subplot(2,1,1);
plot([1:T]', X, 'b-', CP, zeros(size(CP)), 'rx');
grid;

% Now we have some data in X and it's time to perform inference.
% First, setup the matrix that will hold our beliefs about the current
% run lengths.  We'll initialize it all to zero at first.  Obviously
% we're assuming here that we know how long we're going to do the
% inference.  You can imagine other data structures that don't make that
% assumption (e.g. linked lists).  We're doing this because it's easy.
R = zeros([T+1 T]);

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
maxSize = 1;

P = 1;
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
  H = hazard_func([1:t]');
  
  % Evaluate the growth probabilities - shift the probabilities down and to
  % the right, scaled by the hazard function and the predictive
  % probabilities.
  R(2:t+1,t+1) = R(1:t,t) .* predprobs .* (1-H);
  
  % Evaluate the probability that there *was* a changepoint and we're
  % accumulating the mass back down at r = 0.
  R(1,t+1) = sum( R(1:t,t) .* predprobs .* H );
  
  % Renormalize the run length probabilities for improved numerical
  % stability.
  R(:,t+1) = R(:,t+1) ./ sum(R(:,t+1));

  P = [P; R(t+1,t+1)];
  % Update the parameter sets for each possible run length.
  muT0    = [ mu0    ; (kappaT.*muT + X(t)) ./ (kappaT+1) ];
  kappaT0 = [ kappa0 ; kappaT + 1 ];
  alphaT0 = [ alpha0 ; alphaT + 0.5 ];
  betaT0  = [ beta0  ; betaT + (kappaT .*(X(t)-muT).^2)./(2*(kappaT+1)) ];
  muT     = muT0;
  kappaT  = kappaT0;
  alphaT  = alphaT0;
  betaT   = betaT0;
  
  % Store the maximum, to plot later.
  maxes(t) = find(R(:,t)==max(R(:,t)));
  
  if t > 1 && maxes(t,1) - maxes(t-1,1) < 1
  
      importantPoint = [importantPoint; t];
      
  end
  dim = size(R);
  maxSize = max(maxSize, dim(1)*dim(2));
end

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


iFramesInfo.P = P;
iFramesInfo.R = R;
iFramesInfo.maxes = maxes;
iFramesInfo.importantPoint = importantPoint;
iFramesInfo.maxSize = maxSize;

end