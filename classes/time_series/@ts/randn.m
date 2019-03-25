% RANDN Normally distributed pseudorandom numbers.
%    R = RANDN(N) returns an N-by-N matrix containing pseudorandom values drawn
%    from the standard normal distribution.  RANDN(M,N) or RANDN([M,N]) returns
%    an M-by-N matrix. RANDN(M,N,P,...) or RANDN([M,N,P,...]) returns an
%    M-by-N-by-P-by-... array. RANDN returns a scalar.  RANDN(SIZE(A)) returns
%    an array the same size as A.
% 
%    Note: The size inputs M, N, P, ... should be nonnegative integers.
%    Negative integers are treated as 0.
% 
%    R = RANDN(..., CLASSNAME) returns an array of normal values of the 
%    specified class. CLASSNAME can be 'double' or 'single'.
% 
%    R = RANDN(..., 'like', Y) returns an array of normal values of the
%    same class as Y.
% 
%    The sequence of numbers produced by RANDN is determined by the settings of
%    the uniform random number generator that underlies RAND, RANDN, and RANDI.
%    RANDN uses one or more uniform random values to create each normal random
%    value.  Control that shared random number generator using RNG.
%  
%    Examples:
% 
%       Example 1: Generate values from a normal distribution with mean 1
%        and standard deviation 2.
%          r = 1 + 2.*randn(100,1);
% 
%       Example 2: Generate values from a bivariate normal distribution with
%       specified mean vector and covariance matrix.
%          mu = [1 2];
%          Sigma = [1 .5; .5 2]; R = chol(Sigma);
%          z = repmat(mu,100,1) + randn(100,2)*R;
% 
%       Example 3: Reset the random number generator used by RAND, RANDI, and
%       RANDN to its default startup settings, so that RANDN produces the same
%       random numbers as if you restarted MATLAB.
%          rng('default');
%          randn(1,5)
% 
%       Example 4: Save the settings for the random number generator used by
%       RAND, RANDI, and RANDN, generate 5 values from RANDN, restore the
%       settings, and repeat those values.
%          s = rng
%          z1 = randn(1,5)
%          rng(s);
%          z2 = randn(1,5) % z2 contains exactly the same values as z1
% 
%       Example 5: Reinitialize the random number generator used by RAND,
%       RANDI, and RANDN with a seed based on the current time.  RANDN will
%       return different values each time you do this.  NOTE: It is usually
%       not necessary to do this more than once per MATLAB session.
%          rng('shuffle');
%          randn(1,5)
% 
%    See <a href="matlab:helpview([docroot '\techdoc\math\math.map'],'update_random_number_generator')">Replace Discouraged Syntaxes of rand and randn</a> to use RNG to replace
%    RANDN with the 'seed' or 'state' inputs.
% 
%    See also RAND, RANDI, RNG, RANDSTREAM, RANDSTREAM/RANDN
%
%    Reference page in Doc Center
%       doc randn
%
%    Other functions named randn
%
%       codistributed/randn        distributed/randn    RandStream/randn
%       codistributor1d/randn      gpuArray/randn       ts/randn
%       codistributor2dbc/randn
%