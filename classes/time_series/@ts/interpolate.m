%--- help for sde/interpolate ---
%
% INTERPOLATE Stochastic interpolation (Brownian bridge)
% 
%  Syntax:
% 
%    [XT, T] = interpolate(MDL, T, Paths)
%    [XT, T] = interpolate(MDL, T, Paths, 'Name1', Value1, ...)
% 
%  Description:
% 
%    Perform a Brownian interpolation into a user-specified time series array
%    based on a piecewise constant Euler sampling approach.
% 
%    Consider a vector-valued stochastic differential equation (SDE) of the 
%    form
% 
%      dX(t) = F(t,X(t))dt + G(t,X(t))dW(t)
% 
%    where X(t) is an NVARS x 1 state vector, F(t,X(t)) is an NVARS x 1 drift 
%    rate vector-valued function, G(t,X(t)) is an NVARS x NBROWNS diffusion 
%    rate matrix-valued function, and dW(t) is an NBROWNS x 1 Brownian motion 
%    vector. 
% 
%    Given a user-specified time series array associated with the above SDE, 
%    this function performs a stochastic interpolation by sampling from a
%    conditional Gaussian distribution. This sampling technique is sometimes 
%    referred to as a Brownian bridge.
% 
%  Inputs:
% 
%    MDL - Stochastic differential equation model created with the SDE 
%      constructor. See SDE for details.
% 
%    T - NTIMES element vector of interpolation times. The length of T 
%      determine the number of rows in the interpolated output time series 
%      (see XT below).
% 
%    Paths - NPERIODS x NVARS x NTRIALS 3-D time series array of sample paths 
%      of correlated state variables. Each row of Paths is the transpose of 
%      the state vector "X(t)" at time "t" for a given trial. Paths is the 
%      initial time series array into which the Brownian interpolation is 
%      performed.
% 
%  Optional Inputs: 'Name1', Value1, ... is a variable length list of parameter 
%      name/value pairs. The parameter name is specified as a character string, 
%      followed by the corresponding parameter value. Parameter name/value pairs 
%      may be specified in any order; names are case-insensitive and partial 
%      string matches are allowed provided no ambiguities exist. Valid parameter 
%      names are as follows:
% 
%    Times - Vector of monotonically increasing observation times associated 
%      with the time series input Paths. The default is a zero-based, unit-
%      increment column vector of length NPERIODS. 
% 
%    Refine - Scalar logical flag indicating whether the interpolation times 
%      requested by the user (see T above) are used to refine the interpolation 
%      as new information becomes available. If Refine is FALSE, then the 
%      interpolation is based only on the state information specified in Paths. 
%      If Refine is TRUE, then all new interpolated states are inserted into 
%      the existing Paths array as they become available, thereby refining the 
%      interpolation grid available at subsequent interpolation times for the 
%      duration of the current trial. The default is FALSE (no refinement).
% 
%    Processes - Function, or cell array of functions, indicating a sequence 
%      of end-of-period processes or state vector adjustments of the form 
%      X(t) = P(t,X(t)). Processing functions are applied at the end of each 
%      interpolation time, and must accept the current interpolation time "t" 
%      followed by the current state vector "X(t)", and return a state vector 
%      that may be an adjustment to the input state. If more than one processing
%      function is specified, the functions are invoked in the order found in 
%      the cell array. By default no processing or adjustments are made.
% 
%  Outputs:
% 
%    XT - NTIMES x NVARS x NTRIALS 3-D time series array of interpolated state
%      variables formed by interpolating into the input Paths time series array. 
%      Each row of XT is the transpose of the interpolated state vector "X(T)" 
%      at time "T" for a given trial.
% 
%    T - NTIMES x 1 column vector of interpolation times associated with the 
%      output time series (see XT above). If the input interpolation time 
%      vector T contains no missing observations (NaNs), this output is the 
%      same time vector as the input T described above; however, any NaNs are 
%      removed, thereby reducing the length of T and the number of rows of XT.
% 
%  Notes:
% 
%    o All model parameters are assumed piecewise constant, evaluated from 
%      the most recent observation time in Times that precedes a specified
%      interpolation time in T. This is consistent with the Euler approach 
%      of Monte Carlo simulation. 
% 
%    o In the event an interpolation time falls outside the interval spanned 
%      by Times, the time series is extrapolated by an Euler simulation using
%      the nearest available observation. 
% 
%    o The user-defined time series Paths, and corresponding observation Times,
%      must be fully observed (no missing observations indicated by NaNs). 
% 
%    o This method assumes that the user-specified time series array Paths is
%      associated with the SDE object (e.g., the Times/Paths input pair is
%      the result of an initial course-grained simulation). However, the 
%      interpolation ignores the initial conditions of the SDE object (i.e., 
%      MDL.StartTime and MDL.StartState), allowing the Times/Paths input series 
%      specified by the user to take precedence.
% 
%  See also SDE/SIMULATE, SDE.
%
%    Other functions named interpolate
%
%       ts/interpolate
%