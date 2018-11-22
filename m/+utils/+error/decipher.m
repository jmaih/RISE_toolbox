function msgout=decipher(code)
% decipher - interprets error codes
%
% ::
%
%
%   msgout=decipher(code)
%
% Args:
%
%    code (scalar | vector): scalar or vector of error codes returned by
%    RISE.
%
% Returns:
%    :
%
%    - **msgout** [char|cellstr]: explanation of the return code as a char if
%    the input is a scalar or as a cellstr if input is a vector.
%
% Note:
%
% Example:
%
%    See also:

n=numel(code);

if n>1
    
    msgout=cell(n,1);
    
    for ii=1:n
        
        msgout{ii}=decipher(code(ii));
        
    end
    
    return
    
end

msg='';

switch code
    % ====== evaluating the system ====== %
    case 0
    
    case 1
        
        msg='Steady state could not solve';
    
    case 11
        
        msg='Complex or nan steady state residuals';
    
    case 2
        
        msg='Nans in Jacobian';
    
    case 3
        
        msg='Problem in transition matrix';
    
    case 4
        
        msg='Parameter restrictions violated';
    
    case 5
        
        msg='definitions are nan or inf or imaginary';
    
    case 6
        
        msg='Nans or Infs in planner objective';
    
    case 7
        
        msg='bounds or restrictions violations';
        
        % ====== solving the system ====== %
    case 21
        
        msg='Maximum Number of iterations reached or multiple solutions';
    
    case 22
        
        msg='Nans in solution or no solution';
    
    case 223
        
        msg='Complex solution';
    
    case 23
        
        msg='Explosion limit reached';
    
    case 24
        
        msg='Explosive solution';
    
    case 25
        
        msg='System unstable';
    
    case 26
        
        msg='VAR approximation to the DSGE failed';
    
    case 261
        
        msg='The DSGE prior is not proper';
    
    case 27
        
        msg='Cannot take a log-expansion of a variable whose steady state is close to 0';
        
    case 28
        
        msg='For this solver the transition matrix must be diagonal';
        
        % ====== filtering and likelihood evaluation ====== %
    case 301
        
        msg='Maximum number of iterations reached in Lyapunov solution';
    
    case 302
        
        msg='Nans in Lyapunov solution';
    
    case 3002
        
        msg='Nans or Inf in initial conditions for the state vector';
    
    case 30002
        
        msg='Nans or Inf in initial covariance matrix';
    
    case 300002
        
        msg='Nans or Inf or complex likelihood';
    
    case 303
        
        msg='Explosion limit reached in Lyapunov solution';
    
    case 304
        
        msg='Failure of covariance in the determination of sigma points';
    
    case 305
        
        msg='Covariance of forecast errors not positive definite';
    
    case 306
        
        msg='unlikely parameter vector';
    
    case 307
       
        msg='Nan, Inf or complex in log prior';
    
    case 308
        
        msg='Inconsistent ergodic probabilities (nans or sum different from 1)';
    
    case 309
        
        msg='error in the computation of user-defined endogenous priors';
        
        % ====== optimization ====== %
    case 401
        
        msg='optimization failed';
    
    case 402
        
        msg='finding hyperparameters failed';
        
        % ====== loading the data ====== %
    case 500
        
        msg='no actual data or simulation provided for filtering/estimation';
        
        % ====== solving other systems ====== %
    case 201
        
        msg='Non-DSGE: Maximum Number of iterations reached or multiple solutions';
    
    case 202
        
        msg='Non-DSGE: Nans in solution or no solution';
        
        % ====== writing functions to disc ====== %
    case 603
        
        msg='code2file: empty cell';
        
    case 604
        
        msg='code2file: content of cell not a function handle';
        
        % ====== Forecasting ====== %
    case 701
        
        msg='forecast: constraint violation';
        
    case 702
        
        msg='rank deficiency in null and column spaces';
        
    case 703
        
        msg='Simulation: No feasible path';
        
    case 704
        
        msg='Simulation: invalid simulation (complex solutions OR need pruning?)';
        
        % ====== VAR ====== %
    case 801
        
        msg='No suitable rotation could be found';
        
    otherwise
        
        error([mfilename,':: Unknown error code'])
        
end

if nargout
    
    msgout=msg;
    
elseif ~isempty(msg)
    
    disp(msg)
    
end