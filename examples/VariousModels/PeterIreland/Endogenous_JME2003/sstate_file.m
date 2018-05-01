function [y,newp,retcode]=sstate_file(obj,y,p,d,id) %#ok<INUSL>
% sstate_model -- shows the way of writing a RISE steady state file
%
% ::
%
%
%   [y,newp,retcode]=sstate_model(obj,y,p,d,id)
%
% Args:
%
%    - **obj** [rise|dsge]: model object (not always needed)
%
%    - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
%    - **p** [struct]: parameter structure
%
%    - **d** [struct]: definitions
%
%    - **id** [vector]: location of the variables to calculate
%
% Returns:
%    :
%
%      CASE 1: one input argument
%
%    - **y** [cell array]: list of the variables for which the steady state
%    will be calculated within the steady state function
%
%    - **newp** [cell array]: List of the parameters calculated inside the
%    steady state function
%
%    - **retcode** [struct]: possible fields are "imposed", "unique", "loop".
%    The default value for all of those is false.
%      - "imposed": This tells RISE not to check that this is actually solves
%          the steady state. Hence, RISE will attempt to approximate the model
%          around the chosen point
%      - "unique": This tells RISE that the steady state is the same across
%          all regimes. RISE will call the function only once but instead of
%          using just any parameter vector, it will use the ergodic mean of
%          the parameter vector (mean across all regimes).
%      - "loop": This tells RISE that if the user was given the steady state
%          values for some of the variables, he would be able to compute the
%          steady state for the remaining variables. RISE will then exploit
%          this information to optimize over the variables that the user needs
%          for computing the steady state.
%
%      CASE 2: More than one input argument
%
%    - **y** []: endo_nbr x 1 vector of updated steady state
%
%    - **newp** [struct]: structure containing updated parameters if any
%
%    - **retcode** [0|number]: return 0 if there are no problems, else return
%      any number different from 0
%
% Note:
%
%    - If the user knows the steady state, it is always an advantage. If the
%    steady state is computed numerically, we don't know whether it is unique
%    or not. Not that it really matters but... some economists have a strong
%    aversion towards models with multiple equilibria.
%
%    - If the user does not know the solution for all the variables in the
%    steady state, it is a good idea to take a log-linear approximation for
%    the variables that potentially have a nonzero steady state. Hence the
%    user should give that information to RISE.
%
%    - One can potentially improve on the above point by explicit bounds on
%    the variables. But this is not implemented.
%
%    - An alternative that potentially avoids taking a loglinearization is to
%    to reset the values proposed by the optimizer whenever they are in a bad
%    region. It is unclear whether this always works.
%
%    - So be on the safe side, i.e. don't do like me: compute your steady
%    state analytically.
%
% Example:
%
%    See also:

retcode=0;
if nargin==1
    
    % list the variables for which we compute the steady state. Not all
    % variables need to be listed. The variables NOT in this list
    % automatically get a value of 0. RISE is responsible for calculating
    % the steady state of any variable it may have created. The user does
    % not have to worry about those.
    %----------------------------------------------------------------------
    y={'X','A','E','Z','V','MU','PAI','R','Q','LAMBDA','XI','C','M','Y',...
        'K','I','H','W','D','N','LC','LI','LM','LPI','LR'};

	% list of parameters herein calculated
	%-------------------------------------
	newp={};
	
    % flags on the calculation
    %--------------------------
    retcode=struct('unique',true,'imposed',true);
else
    % no parameters to update
    %-------------------------
    newp=[];
    
    % create definitions (auxiliary parameters) that are useful for solving
    % the steady state
    %----------------------------------------------------------------------
    z_ss=p.z_ss_trans*10000;
% 	phi_p = 100*abs(p.phi_p_trans);
% 	phi_k = 100*abs(p.phi_k_trans);
	mu_ss = 1 + p.mu_ss_trans;

    X=1;
    A=1;
    E=p.e_ss;
    Z=z_ss;
    V=1;
    % mu determined by policy
    MU=mu_ss;
    PAI=MU/p.g;
    R=MU/p.beta;
    Q=p.g/p.beta-1+p.delta;
    % lambda comes here
    l1ss=1/(p.theta/(p.theta-1)-(p.g-1+p.delta)*p.alpha/Q);
    l2ss=1/(1+E*(R/(R-1))^(p.gam-1));
    l3ss=(1-p.alpha)*Z*((p.theta-1)/p.theta)^(1/(1-p.alpha))*(p.alpha/Q)^(p.alpha/(1-p.alpha));
    LAMBDA=(p.eta+(1-p.alpha)*l1ss*l2ss)/l3ss;
    XI = (p.theta-1)/p.theta*LAMBDA;
    C = 1/(1+E*(R/(R-1))^(p.gam-1))*(1/LAMBDA);
    M = E*(R/(R-1))^p.gam*C;
    Y = 1/(1-(p.g-1+p.delta)*(p.theta-1)/p.theta*p.alpha/Q)*C;
    K = (p.theta-1)/p.theta*p.alpha*Y/Q;
    I = (p.g-1+p.delta)*K;
    H=1/Z*(Y/K^p.alpha)^(1/(1-p.alpha));
    W=(1-p.alpha)*(p.theta-1)/p.theta*Y/H;
    D = Y-W*H-Q*K;
    N=Y/H;
    LC=log(C);
    LI=log(I);
    LM=log(M);
    LPI=log(PAI);
    LR=log(R);
    ys=[X A E Z V MU PAI R Q LAMBDA XI C M Y K I H W D N LC LI LM LPI LR].';
    
    % check the validity of the calculations
    %----------------------------------------
    if ~utils.error.valid(ys)||any(abs(ys)>100000)
        retcode=1;
    else
        % push the calculations
        %----------------------
        y(id)=ys;
    end
end

end