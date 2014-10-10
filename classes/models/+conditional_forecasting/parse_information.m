function conditions=parse_information(conditional_info,endo_names,...
    exo_names,forecast_start_date)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(conditional_info)
    conditions=[];
    return
end
endo_vars={};
exo_vars={};
min_date=inf;
max_date=-inf;
fields=fieldnames(conditional_info);
main_fields={'lower_bound','upper_bound','central_tendency'};

get_names_and_dates();

date_span=min_date:max_date;
nobs=numel(date_span);
all_vars=[endo_vars,exo_vars];

conditions=rebuild_conditions();

format_conditions();

    function format_conditions()
        % holes are allowed in LB and UB only in places where the CT has
        % holes
        nan_ct=isnan(conditions.central_tendency);
        conditions.lower_bound(isnan(conditions.lower_bound) & ~nan_ct)=-inf;
        conditions.upper_bound(isnan(conditions.upper_bound) & ~nan_ct)=inf;
        % now format the conditions
        conditions.OMG_endo=[];
        conditions.OMG_exo=[]; % for the moment, I do not condition on a particular covariance matrix
        nendo=numel(endo_vars);
        nexo=numel(exo_vars);
        conditions.restrictions_id=nan(1,nendo+nexo);
        conditions.restrictions_id(1:nendo)=locate_variables(endo_vars,endo_names);
        conditions.restrictions_id(nendo+(1:nexo))=locate_variables(exo_vars,exo_names);
        conditions.endogenous_vars=false(1,nendo+nexo);
        conditions.endogenous_vars(1:nendo)=true;
        conditions.varnames=all_vars;
    end
    function conditions=rebuild_conditions()
        tmp=nan(nobs,numel(all_vars));
        conditions=struct('lower_bound',tmp,'upper_bound',tmp,'central_tendency',tmp);
        for ii=1:numel(fields)
            new_ts=conditional_info.(fields{ii});
            dn=new_ts.date_numbers;
            date_locs=find(date_span==dn(1));
            locs=locate_variables(new_ts.varnames,all_vars);
            conditions.(fields{ii})(date_locs-1+(1:numel(dn)),locs)=double(new_ts);
        end
        % check start date of the conditions with respect to the start of
        % the forecasts
        if min_date<forecast_start_date
            error([mfilename,':: conditional information cannot start before the beginning of the forecast period'])
        end
        % add the missing conditioning information dates...
        addendum=numel(forecast_start_date:min_date)-1;
        addendum=nan(addendum,numel(all_vars));
        for ii=1:numel(main_fields)
            conditions.(main_fields{ii})=[addendum;conditions.(main_fields{ii})];
        end
    end
    function get_names_and_dates()
        for ii=1:numel(fields)
            if ismember(fields{ii},main_fields)
                new_ts=ts.collect(conditional_info.(fields{ii}));
                if new_ts.NumberOfPages~=1
                    error([mfilename,':: conditional database ',fields{ii},' must have one page only'])
                end
                if ii==1
                    frequency=new_ts.frequency;
                else
                    if ~strcmp(frequency,new_ts.frequency)
                        error('data should have the same frequency')
                    end
                end
                min_date=min(min_date,new_ts.date_numbers(1));
                max_date=max(max_date,new_ts.date_numbers(end));
                if ~isa(new_ts,'ts')
                    error([mfilename,'field ',fields{ii},' must be a ts object'])
                end
                if isempty(new_ts.varnames{1})
                    error([mfilename,':: conditional information variables should have names'])
                end
                for oo=1:new_ts.NumberOfVariables
                    vo=new_ts.varnames{oo};
                    if ismember(vo,endo_names)
                        endo_vars=union(endo_vars,vo);
                    elseif ismember(vo,exo_names)
                        exo_vars=union(exo_vars,vo);
                    else
                        error([mfilename,':: ',vo,' not recognized as an endogenous or an exogenous variable'])
                    end
                end
            else
                error([mfilename,':: conditional info fields must elements of lower_bound, upper_bound, central_tendency'])
            end
        end
    end
end
