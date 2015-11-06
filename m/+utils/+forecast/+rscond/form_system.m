function [RS,yhats]=form_system(M,ycond,econd,nsteps,options)

%   - **forecast_conditional_hypothesis**
%   [{'default'}|'ncp'/'econ'|'nas'/'econ']: in dsge models in which agents
%   have information beyond the current period, this option determines the
%   number of periods of shocks need to match the restrictions:
%       - Hypothesis **'default'** assumes that irrespective of how
%           many periods of conditioning information are remaining, agents
%           always receive information on the same number of shocks.
%       - Hypothesis **ncp or econ** assumes there are as many shocks periods as
%           the number of the number of conditioning periods
%       - Hypothesis **nas** assumes there are as many shocks periods as
%           the number of anticipated steps. It has now been merged with
%           econ
%

% forecast_conditional_hypothesis econ seems to create some inaccuracies
% - Different inversions of the matrices give different solutions despite
% the fact that the matrices are full ranked. 
% - This also leads to unrealistically-sized shocks... 
% - We have checked that the restrictions remain well enforced albeit with
% negligeable inaccuracies
% - What I have not thoroughly checked is whether the other inverses lead
% could lead to smaller shocks. The first investigations reveal that this
% is not the case.
% - The problem does not depend on whether we have hard or soft
% restrictions, and it is probably worsened in the latter case.

hypotheses={'default','jma','ncp','econ','nas'};
defaults={
    'debug',false,@(x)islogical(x),'debug must be a logical'
    
    'forecast_conditional_hypothesis','default',@(x)ismember(x,hypotheses),...
    ['forecast_conditional_hypothesis must be one of ',cell2mat(strcat(hypotheses,'|'))]
    
    'nsteps',[],@(x)isnumeric(x)&&isscalar(x)&&x>0 && floor(x)==ceil(x),...
    'nsteps must be a positive integer'
    };

if nargin==0
    RS=cell2struct(defaults(:,2),defaults(:,1),1);
    yhats=defaults;
    return
else
    narginchk(4,5)
    if nargin<5
        options=[];
    end
end

if isempty(options)
    options=struct();
end

options=parse_arguments(defaults,options);
% stack R and S
%-----------------
RS=[M.R;M.S];
rs_cols=size(RS,2);

% Trim the data to have a length of nsteps by chopping or adding nans
%--------------------------------------------------------------------
[conddatay,ry,cy]=allocate(ycond);
[conddatae,~,ce]=allocate(econd,true);

% Demean the conditions in 3 pages: central,lb,ub
%-------------------------------------------------
yhat=bsxfun(@minus,conddatay,reshape(M.ufkst+M.const,ry,nsteps));
yhat=reshape(yhat,[],3);
% vectorize and add the shocks conditions
%----------------------------------------
ehat=reshape(conddatae,[],3);
yhats=[yhat;ehat]; % central,lb,ub

% remove the bad rows: rows with 3 nan
%--------------------------------------
good=~all(isnan(yhats),2);
yhats=yhats(good,:);
if any(isnan(yhats(:)))
    error('nan remaining in conditions: central,lb,ub, for an observations should all be nan or not')
end
RS=RS(good,:);

% remove bad/extra columns
%---------------------------
% This also has to be squared with the randomness of the shocks above
% because in case of random shocks beyond the number of columns of RS
% (chopped or not), they will be taken into account when re-forecasting
% RS =RS(:,1:?)
% At this point, we have length_data=cy=ce=nsteps
length_data=max(cy,ce);
hypo=upper(options.forecast_conditional_hypothesis);
switch hypo
    case {'NCP','ECON','NAS'}
        % This will always match with nsteps in this formulation. It is the
        % responsibility of the user to make sure that the number of
        % observations is equal to the number of steps
        maxcols=length_data*M.nshocks;
%     case {'NAS'}
%         nap=k+1;
%         maxcols=nap*M.nshocks;
%         % This is good for one-step forecast...
%         if ~isequal(length_data,nap)
%             error([mfilename,':: for the NAS assumption, you need # anticipated steps = # conditioning periods'])
%         end
    case {'DEFAULT','JMA'}
        % no chopping
        maxcols=rs_cols;
    otherwise
        error([mfilename,':: Unknown option for the anticipation forecast_conditional_hypothesis'])
end

rs_cols=min(rs_cols,maxcols);

RS=RS(:,1:rs_cols);



    function [A,r,c]=allocate(xcond,isshock)
        if nargin<2
            isshock=false;
        end
        if size(xcond.data,3)==1
            % hard conditions: lower-bound equals upper bound equals mean
            %------------------------------------------------------------
            xcond.data=xcond.data(:,:,ones(1,3));
        end
        [r,c,pages]=size(xcond.data);
        if pages~=3
            error('wrong number of pages')
        end
        if options.debug
            if c>nsteps
                % no problem
                disp('More conditioning periods than the forecast horizon')
            elseif c<nsteps
                % we are conditioning on few periods than the forecast horizon
                % and that is OK
                disp('Fewer conditioning periods than the forecast horizon')
            else
                disp('As many conditioning periods as the forecast horizon')
            end
        end
        if isshock
            A=xcond.data;
            if r>0
                rstar=M.nshocks;
                missing=rs_cols-rstar*c;
                cplus=missing/rstar;
                if missing>0
                    % add nans
                    A=cat(2,A,nan(r,cplus,3));
                elseif missing<0
                    % remove extra columns
                    A=A(:,1:c+cplus,:);
                end
            end
        else
            nhard=min(nsteps,c);
            cutoff=min(nhard,nsteps);
            A=nan(r,nsteps,3);
            A(:,1:cutoff,:)=xcond.data(:,1:cutoff,:);
        end
    end
end