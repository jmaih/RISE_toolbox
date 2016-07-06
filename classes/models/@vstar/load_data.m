function [obj,issue,retcode]=load_data(obj,varargin)
% load_data - loads data for estimation and for forecasting
%
% Syntax
% -------
% ::
%
%   [obj,issue,retcode]=load_data(obj)
%
%   [obj,issue,retcode]=load_data(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: scalar or vector or RISE model objects
%
% - **data** [ts|struct]: data to load
%
% - **data_demean** [true|{false}]: flag for demeaning the data
%
% Outputs
% --------
%
% - **obj** [rise|dsge|rfvar|svar]: scalar or vector or RISE model objects
%
% - **issue** [{''}|char]: message describing any issue related to the
% loading of the data
%
% - **retcode** [{0}|positive integer]: if retcode is different from 0,
% then there is a problem
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:


if isempty(obj)
    
    obj=struct('data',ts,...
        'data_demean',false);
    
    return
    
end

nobj=numel(obj);

issue='';

retcode=0;

for ii=1:nobj
            
        if obj(ii).observables.number(1)==0
            
            error([mfilename,':: List of observables not provided for model ''',obj(ii).filename,''])
            
        end
        
        [obj(ii),issue,rcode]=load_data_intern(obj(ii),varargin{:});
        
        if nobj>1 && ~isempty(issue)
            
            warning([mfilename,':: problem for model ''',obj(ii).filename,''': ',issue]) %#ok<WNTAG>
            
        end
        
        retcode=max(retcode,rcode);
        
    
end

end

function [obj,issue,retcode]=load_data_intern(obj,varargin)

issue='';

retcode=0;

obj=set(obj,varargin{:});

% check wether the data are provided
obj.options.data=ts.collect(obj.options.data);

data_provided=obj.options.data.NumberOfVariables>0;

if ~data_provided
    
    error('No data to fit the model')
    
end

allvars=obj.options.data.varnames;

pages=1;

[verdier,obj.options.estim_start_date,obj.options.estim_end_date]=...
    utils.time_series.data_request(obj.options.data,allvars,...
    obj.options.estim_start_date,obj.options.estim_end_date,pages);

[verdier,obj.options.estim_start_date,obj.options.estim_end_date]=sweep_nans(...
    verdier,obj.options.estim_start_date,obj.options.estim_end_date);

% now account for lags that appear in the transition variables
%-------------------------------------------------------------
second_sweep()

if obj.options.data_demean
    
    verdier=bsxfun(@minus,verdier,mean(verdier,2));
    
end
    
    % add further description fields
    %-------------------------------
    obj.data=verdier;
    
    obj.data_are_loaded=true;

    function second_sweep()
        % rebuild the dataset and set the data according to the description in the
        % observables...
        newdata=ts(obj.options.estim_start_date,verdier.',allvars);
        
        newdata=pages2struct(newdata);
        
        batch=ts.empty;
        
        for ivar=1:sum(obj.observables.number)
            
            vname=obj.observables.name{ivar};
            
            lag=0;
            
            lp=find(vname=='{');
            
            if isempty(lp)
                
                lp=find(vname=='(');
                
            end
            
            if ~isempty(lp)
                
                rp=find(vname=='}');
                
                if isempty(rp)
                    
                    rp=find(vname==')');
                    
                end
                
                if isempty(rp)
                    
                    error('wrong variable name in observables')
                    
                else
                    
                    lag=str2double(vname(lp+1:rp-1));
                    
                    vname=vname(1:lp-1);
                    
                end
                
            end
            
            batch=[batch,newdata.(vname){lag}];
            
        end
        
        verdier=double(batch).';
        
        [verdier,obj.options.estim_start_date,obj.options.estim_end_date]=...
            sweep_nans(verdier,batch.start);
    end

end

function [verdier,start_date,end_date]=sweep_nans(verdier,start_date,~)
% trim for missing/nans
bad=any(isnan(verdier),1);

first_good=find(bad==0,1,'first');

last_good=find(bad==0,1,'last');

verdier=verdier(:,first_good:last_good);

if isempty(verdier)
    
    error('Insufficient data for estimation')
    
end

% re-adjust the dates
%--------------------
start_date=obs2date(start_date,first_good);

end_date=obs2date(start_date,size(verdier,2));

end