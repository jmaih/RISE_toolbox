function [M,ufkst,states,PAI,TT,Q]=stochastic_impact_cumulator(model,y0,nsteps,...
y_pos,e_pos,states)
% STOCHASTIC_IMPACT_CUMULATOR -- creates impact matrix for contemporaneous
% and future shocks
%
% ::
%
%
%   M=STOCHASTIC_IMPACT_CUMULATOR(model,y0,nsteps)
%
%   M=STOCHASTIC_IMPACT_CUMULATOR(model,y0,nsteps,y_pos)
%
%   M=STOCHASTIC_IMPACT_CUMULATOR(model,y0,nsteps,y_pos,e_pos)
%
%   M=STOCHASTIC_IMPACT_CUMULATOR(model,y0,nsteps,y_pos,e_pos,states)
%
%   [M,ufkst,states,PAI,TT,Q]=STOCHASTIC_IMPACT_CUMULATOR(...)
%
% Args:
%
%    - **model** [struct]:
%      - **T** [1 x h cell]: solution of the model,
%      T{regime}=[Ty,Tsig,Te_0,...Te_k]
%      - **sstate** [1 x h cell]: steady state in each regime
%      - **state_cols** [vector|{1:ny}]: location of state endogenous
%      variables in y0 (see below).
%      - **k** [scalar]: anticipation horizon (Beyond the current period)
%      - **Qfunc** [empty|function_handle]: endogenous transition matrix
%
%    - **y0** [ny x 1 vector]: initial conditions
%
%    - **nsteps** [integer]: number of forecast steps
%
%    - **y_pos** [vector]: location of restricted endogenous variables.
%
%    - **e_pos** [vector]: location of restricted shocks
%
%    - **states** [empty|nsteps x 1 vector]: states visited for each forecast
%    step.
%
% Returns:
%    :
%
%    - **M** [struct]:
%      - **R** [matrix]: convoluted restrictions of shocks stemming from the
%      restrictions on endogenous variables
%      - **ufkst** [matrix]: unconditional forecasts for the restricted
%      endogenous variables (Excluding the initial conditions!!!).
%      - **const** [vector]: impact of the constant (steady state + risk) for
%      the restricted endogenous variables.
%      - **S** [matrix]: direct restrictions on shocks
%      - **nshocks** [integer]: number of shocks
%      - **ny** [integer]: number of endogenous variables
%
%    - **ufkst** [ny x (nsteps+1) matrix]: Unconditional forecasts mean (with
%    initial condition at the beginning!!!)
%
%    - **states** [nsteps x 1 vector]: regimes visited for each forecast step.
%
%    - **PAI** [h x nsteps matrix]: choice probabilities of regimes for each
%    step
%
%    - **TT** [matrix]: convolution of autoregressive terms for the
%    restrictions
%
%    - **Q** [h x h x nsteps array]: Time series of transition matrices
%
% Note:
%
% Example:
%
%    See also:

% - Here we don't need to know anything about how long we have data for.
% After building the matrix, we can use the hypotheses to chop off the
% extra columns
% - admissibility, then, will reflect the fact that the maximum horizon of
% shocks may be constrained by the availability of the data, depending on
% the hypothesis entertained.
% - It will also be possible to kill the non-active shocks although,
% perhaps, that could already be baked into the solution...
% - the steady state can no longer be dealt with only at the end because it
% may change from one state to another
% - What is the most likely combination of shocks (and states) that
% minimizes some variance... This could be an optimization problem...
% whereby we optimize over both states and shocks. But then, the problem
% does not have to be linear or conditionally linear. This is a more
% general problem in the sense that we can even add restrictions such as
% zlb and so on.
% - entropy forecasting

% FileInfo = dir('bvar_rise.m')
% http://www.mathworks.com/matlabcentral/answers/33220-how-do-i-get-the-file-date-and-time
% fs2000sims.mat
% altmany at gmail.com

if nargin<6
    
    states=[];
    
    if nargin<5
        
        e_pos=[];
        
        if nargin<4
            
            y_pos=[];
            
        end
        
    end
    
end

if isempty(states)
    
    states=nan(nsteps,1);
    
end

if numel(states)~=nsteps
    
    error('number of states does not match number of steps')
    
end

T=model.T;

[ny,nz]=size(T{1});

sstate=model.sstate;

state_cols=model.state_cols;

if isempty(state_cols)
    
    state_cols=1:ny;
    
end

k=model.k;

Qfunc=model.Qfunc;

h=size(T,2);

nx=numel(state_cols);

nshocks=(nz-(nx+1))/(k+1);

[Ty,Te,~,C]=utils.forecast.rscond.separate_terms(T,sstate,state_cols,k,nshocks);

% const/y{-1}/shk(0)/shk(1)/.../shk(k)/shk(k+1)/.../shk(k+nsteps-1)
ProtoR=cell(1,1+1+k+nsteps);

ProtoR{1}=zeros(ny,1);

ProtoR{2}=y0;

ncols=1+1+(k+nsteps)*nshocks;

nconds=numel(y_pos);

R=zeros(nconds*nsteps,ncols);

is_convolute=nargout>4;

if is_convolute
    
    TT=zeros(nconds*nsteps,ny);
    
    Tconv=eye(ny);
    
end

ufkst=y0(:,ones(1,nsteps+1));

PAI=zeros(h,nsteps);

Q=nan(h,h,nsteps);

for jstep=1:nsteps
    
    % pick the state
    %----------------
    [st,Q(:,:,jstep)]=pick_a_regime();
    
    % deal with the constant
    %------------------------
    ProtoR{1}=C{st}+Ty{st}*ProtoR{1};
    
    % origin
    %--------
    ProtoR{2}=Ty{st}*ProtoR{2};
    
    % convolution of autoregressive terms
    %-------------------------------------
    if is_convolute
        
        convolute()
        
    end
    
    % shocks
    %--------
    iter=2;
    
    for icol=1:k+jstep
        
        iter=iter+1;
        
        if isempty(ProtoR{iter})
            
            ProtoR{iter}=zeros(ny,nshocks);
            
        else
            
            ProtoR{iter}=Ty{st}*ProtoR{iter};
            
        end
        
        which_shock=icol-jstep+1;
        
        if which_shock>0
            
            ProtoR{iter}=ProtoR{iter}+Te{st}{which_shock};
            
        end
        
    end
    
    % process all cells up to iter
    %------------------------------
    tmp=cell2mat(ProtoR(:,1:iter));
    
    R((jstep-1)*nconds+1:jstep*nconds,1:size(tmp,2))=tmp(y_pos,:);
    
    ufkst(:,jstep+1)=tmp(:,2); % unconditional forecasts
    
    % next period's origin
    %----------------------
    y0=ProtoR{1}+ProtoR{2};
    
end

% create rows for the shocks
%-----------------------------
S=do_shocks_conditions();

% format output
%---------------
R=sparse(R);

M=struct('R',R(:,3:end),'ufkst',R(:,2),'const',R(:,1),'S',S,...
    'nshocks',nshocks,'ny',ny);

    function convolute()
        
        Tconv=Ty{st}*Tconv;
        
        TT((jstep-1)*nconds+1:jstep*nconds,:)=Tconv(y_pos,:);
        
    end

    function [reg,Q]=pick_a_regime()
        
        % compute the transition matrix irrespective of whether the regimes
        % are preset or not
        Q=Qfunc(y0);
        
        if isnan(states(jstep))
            
            if h==1
                
                states(jstep)=1;
                
                PAI(:,jstep)=1;
                
            else
                
                if jstep==1
                    
                    [PAI0,retcode]=initial_markov_distribution(Q,true);
                    
                    if retcode
                        
                        warning('ergodic distribution failed...')
                        
                        PAI0=initial_markov_distribution(Q,false);
                        
                    end
                    
                    PAI(:,jstep)=PAI0(:);
                    
                else
                    
                    % draw conditional on yesterday's state
                    PAI(:,jstep)=Q(states(jstep-1),:).';
                    
                end
                
                cp=cumsum(PAI(:,jstep).');
                
                cp=[0,cp];
                
                states(jstep)=find(cp>rand,1,'first')-1;
                
            end
            
        end
        
        reg=states(jstep);
        
    end

    function S=do_shocks_conditions()
        
        S=speye(nshocks*(k+nsteps));
        
        e_good=false(nshocks,1);
        
        e_good(e_pos)=true;
        
        e_good=e_good(:,ones(1,k+nsteps));
        
        S=S(e_good(:),:);
        
    end

end