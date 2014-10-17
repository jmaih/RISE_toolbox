function [loglik,Incr,retcode,Filters]=msre_linear_filter(syst,data_info,data_trend,state_trend,SS,risk,options)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% all rows of Q should sum to 1
% 0: no filters,
% 1: filtering,
% 2: filtering+updating,
% 3: filtering+updating+smoothing
% this decides the initial condition for both the state and the markov chain distributions in the
% the kalman filter.
% now add the defaults for the initialization process
if nargin==0
    if nargout>1
        error([mfilename,':: with no input argument, the number of output arguments cannot exceed 1'])
    end
    filt_options=struct('kf_algorithm','lwz',...%     alternative is kn (Kim and Nelson)
        'kf_tol',1e-20,...
        'kf_filtering_level',3,...
        'kf_riccati_tol',1e-6,...
        'kf_nsteps',1);
    init_options=kalman_initialization();
    defaults=utils.miscellaneous.mergestructures(filt_options,init_options);
    loglik=defaults;
    return
end


% apply the presample to the data
%--------------------------------
data_info.include_in_likelihood(1:options.kf_presample)=false;

% Trim the system matrices if possible
%---------------------------------------
minimum_state_for_estimation()

[init,retcode]=kalman_initialization(syst.T,syst.R,SS,risk,syst.Qfunc,options);

if retcode
    loglik=[];
    Incr=[];
    Filters=[];
else
    [loglik,Incr,retcode,Filters]=msre_kalman_cell(syst,data_info,data_trend,state_trend,init,options);
end
if isempty(loglik)
    loglik=nan;
end
% agenda:
% add back the Kim-Nelson filter
% do the initialization for non-stationary variables
% add the Kim-Nelson filter with k-states carried
% The real-time filter will need to be added as well at a later stage.
% Create the environment for the new deal
% branch the derivatives
% solve the second order

    function minimum_state_for_estimation()
        % this function reduces the state size to accelerate estimation
        % it returns the indices of the union of state variables and observable
        % ones
    % N.B: the endogenous probability script is written as a function of
    % the grand state (i.e. all the endogenous variables) and so, it is
    % dangerous to reduce the state in that case. The drawback is that
    % likelihood computation becomes slower. Amending this would amount to
    % having a separate script for the transition matrix for the specific
    % case where we don't need all the smoothed variables. This will have
    % to be done from the parser, when reading the list of the observables.
    % Or something along those lines...
    
    % IF THE STATE IS SHRUNK, JUST WRITE A WRAPPER THAT WILL INFLATE IT
    % BEFORE IT IS APPLIED TO SYST.QFUNC. IN BRIEF, SYST.QFUNC WILL BE
    % MODIFIED TO SOMETHING LIKE syst.Qfunc=@(x)syst.Qfunc(inflator(x)).
    % The problem is that the inflator should know how to reset the grand
    % vector and the variables entering the transition matrix should be
    % forced to be states...
    return
        if options.kf_filtering_level==0 
            tmp=syst.T;
            h=numel(tmp);
            % The step below is critical for speed, though it may add some noise
            % when computing the filters for all the endogenous variables
            state=any(abs(tmp{1})>1e-9,1);
            for ii=2:h
                state=state | any(abs(tmp{ii})>1e-9,1);
            end
            n=numel(state);
            newstate=state;
            newstate(data_info.varobs_id)=true;
            newobs=false(1,n);
            newobs(data_info.varobs_id)=true;
            newobs(~newstate)=[];
            newobs=find(newobs);
            [~,tags]=sort(data_info.varobs_id);
            newobs(tags)=newobs;
            data_info.varobs_id=newobs;
            for istate=1:h
                SS{istate}=SS{istate}(newstate);
                risk{istate}=risk{istate}(newstate);
                syst.T{istate}=syst.T{istate}(newstate,newstate);
                syst.R{istate}=syst.R{istate}(newstate,:,:);
                if ~isempty(state_trend{1})
                    state_trend{istate}=state_trend{istate}(newstate,:);
                end
            end
        end
    end
end