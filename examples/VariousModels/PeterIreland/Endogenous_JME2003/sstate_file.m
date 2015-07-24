function [y,newp,retcode]=sstate_file(obj,y,p,d,id) %#ok<INUSL>
% sstate_model -- shows the new way of writing a RISE steady state file
%
% Syntax
% -------
% ::
%
%   [y,newp,retcode]=sstate_model(obj,y,p,d,id)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model object (not always needed)
%
% - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
% - **p** [struct]: parameter structure
%
% - **d** [struct]: definitions
%
% - **id** [vector]: location of the variables to calculate
%
% Outputs
% --------
%
% - **y** []: endo_nbr x 1 vector of updated steady state
%
% - **newp** [struct]: structure containing updated parameters if any
%
% - **retcode** [0|number]: return 0 if there are no problems, else return
%   any number different from 0
%
% More About
% ------------
%
% - this is new approach has three main advantages relative to the previous
%   one:
%   - The file is valid whether we have many regimes or not
%   - The user does not need to know what regime is being computed
%   - It is in sync with the steady state model
%
% Examples
% ---------
%
% See also:

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
    % flags on the calculation
    %--------------------------
    newp=struct('unique',true,'imposed',true,'initial_guess',false);
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
    if ~utils.error.valid(ys)
        retcode=1;
    else
        % push the calculations
        %----------------------
        y(id)=ys;
    end
end

end