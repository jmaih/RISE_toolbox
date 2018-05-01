function [obj,issue,retcode]=load_data(obj,varargin)
% load_data - loads data for estimation and for forecasting
%
% ::
%
%
%   [obj,issue,retcode]=load_data(obj)
%
%   [obj,issue,retcode]=load_data(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar]: scalar or vector or RISE model objects
%
%    - **data** [ts|struct]: data to load
%
%    - **data_demean** [true|{false}]: flag for demeaning the data
%
% Returns:
%    :
%
%    - **obj** [rise|dsge|rfvar|svar]: scalar or vector or RISE model objects
%
%    - **issue** [{''}|char]: message describing any issue related to the
%    loading of the data
%
%    - **retcode** [{0}|positive integer]: if retcode is different from 0,
%    then there is a problem
%
% Note:
%
% Example:
%
%    See also:

if isempty(obj)
    
    mydefaults=the_defaults();
    
    if nargout
        
        obj=mydefaults;
        
    else
        
        clear obj
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

nobj=numel(obj);

issue='';

retcode=0;

for ii=1:nobj
    
    if isa(obj(ii),'dsge') && obj(ii).is_optimal_simple_rule_model
        
        if obj(ii).observables.number(1)>0
            
            warning([mfilename,':: data not needed OPTIMAL SIMPLE RULES problems'])
            
        end
        
        obj(ii).data_are_loaded=true;
        
        continue
        
    else
        
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

end

function [obj,issue,retcode]=load_data_intern(obj,varargin)

issue='';

obj=set(obj,varargin{:});

[dataOut,dataset,estim_start_date,estim_end_date,retcode]=...
    data_prerequest(obj,obj.options.data);

% check wether the data are provided
if isempty(dataset)
    
    return
    
end

obj.options.data=dataset;

obj.data=dataOut;
    
obj.data_are_loaded=true;
    
obj.options.estim_start_date=estim_start_date;
    
obj.options.estim_end_date=estim_end_date;

end

function d=the_defaults()

d={'data',ts,@(x)isa(x,'ts')||isstruct(x),...
    'data must be a ts or a structure of ts'
    
    'data_demean',false,@(x)islogical(x),...
    'data_demean must be a logical'
    };

end


