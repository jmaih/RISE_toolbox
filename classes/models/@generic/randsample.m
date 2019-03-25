% RANDSAMPLE Random sample, with or without replacement.
%    Y = RANDSAMPLE(N,K) returns Y as a column vector of K values sampled
%    uniformly at random, without replacement, from the integers 1:N.
% 
%    Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at random,
%    without replacement, from the values in the vector POPULATION.  Y is a
%    vector of the same type as POPULATION.  NOTE:  When POPULATION is a
%    numeric vector containing only non-negative integer values, and it might
%    have length 1, use
% 
%       Y = POPULATION(RANDSAMPLE(LENGTH(POPULATION),K)
% 
%    instead of Y = RANDSAMPLE(POPULATION,K).
% 
%    Y = RANDSAMPLE(N,K,REPLACE) or RANDSAMPLE(POPULATION,K,REPLACE) returns a
%    sample taken with replacement if REPLACE is true, or without replacement
%    if REPLACE is false (the default).
% 
%    Y = RANDSAMPLE(N,K,true,W) or RANDSAMPLE(POPULATION,K,true,W) returns a
%    weighted sample, using positive weights W, taken with replacement.  W is
%    often a vector of probabilities. This function does not support weighted
%    sampling without replacement.
% 
%    Y = RANDSAMPLE(S,...) uses the random number stream S for random number 
%    generation.  RANDSAMPLE uses the MATLAB default random number stream by
%    default.
% 
%    Examples:
% 
%    Draw a single value from the integers 1:10.
%       n = 10;
%       x = randsample(n,1);
% 
%    Draw a single value from the population 1:n, when n > 1.
%       y = randsample(1:n,1);
% 
%    Generate a random sequence of the characters ACGT, with
%    replacement, according to specified probabilities.
%       R = randsample('ACGT',48,true,[0.15 0.35 0.35 0.15])
% 
%    See also RAND, RANDPERM, RANDSTREAM.
%
%    Reference page in Doc Center
%       doc randsample
%
%    Other functions named randsample
%
%       generic/randsample
%