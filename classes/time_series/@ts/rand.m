% RAND Uniformly distributed pseudorandom numbers.
%    R = RAND(N) returns an N-by-N matrix containing pseudorandom values drawn
%    from the standard uniform distribution on the open interval(0,1).  RAND(M,N)
%    or RAND([M,N]) returns an M-by-N matrix.  RAND(M,N,P,...) or
%    RAND([M,N,P,...]) returns an M-by-N-by-P-by-... array.  RAND returns a
%    scalar.  RAND(SIZE(A)) returns an array the same size as A.
% 
%    Note: The size inputs M, N, P, ... should be nonnegative integers.
%    Negative integers are treated as 0.
% 
%    R = RAND(..., CLASSNAME) returns an array of uniform values of the 
%    specified class. CLASSNAME can be 'double' or 'single'.
% 
%    R = RAND(..., 'like', Y) returns an array of uniform values of the 
%    same class as Y.
% 
%    The sequence of numbers produced by RAND is determined by the settings of
%    the uniform random number generator that underlies RAND, RANDI, and RANDN.
%    Control that shared random number generator using RNG.
% 
%    Examples:
% 
%       Example 1: Generate values from the uniform distribution on the
%       interval (a, b).
%          r = a + (b-a).*rand(100,1);
% 
%       Example 2: Use the RANDI function, instead of RAND, to generate
%       integer values from the uniform distribution on the set 1:100.
%          r = randi(100,1,5);
% 
%       Example 3: Reset the random number generator used by RAND, RANDI, and
%       RANDN to its default startup settings, so that RAND produces the same
%       random numbers as if you restarted MATLAB.
%          rng('default')
%          rand(1,5)
% 
%       Example 4: Save the settings for the random number generator used by
%       RAND, RANDI, and RANDN, generate 5 values from RAND, restore the
%       settings, and repeat those values.
%          s = rng
%          u1 = rand(1,5)
%          rng(s);
%          u2 = rand(1,5) % contains exactly the same values as u1
% 
%       Example 5: Reinitialize the random number generator used by RAND,
%       RANDI, and RANDN with a seed based on the current time.  RAND will
%       return different values each time you do this.  NOTE: It is usually
%       not necessary to do this more than once per MATLAB session.
%          rng('shuffle');
%          rand(1,5)
% 
%    See <a href="matlab:helpview([docroot '\techdoc\math\math.map'],'update_random_number_generator')">Replace Discouraged Syntaxes of rand and randn</a> to use RNG to replace
%    RAND with the 'seed', 'state', or 'twister' inputs.
% 
%    See also RANDI, RANDN, RNG, RANDSTREAM, RANDSTREAM/RAND,
%             SPRAND, SPRANDN, RANDPERM.
%
%    Reference page in Doc Center
%       doc rand
%
%    Other functions named rand
%
%       codistributed/rand        distributed/rand    RandStream/rand
%       codistributor1d/rand      gpuArray/rand       ts/rand
%       codistributor2dbc/rand
%