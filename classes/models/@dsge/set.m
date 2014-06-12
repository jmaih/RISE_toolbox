function obj=set(obj,varargin)
if isempty(obj)
    obj=set@rise_generic(obj);
    % add other fields here if necessary
else
    % obj=set(obj,'shock_horizon',struct('shock1',2,'shock3',4))
    % obj=set(obj,'shock_horizon',5)
    shock_horizon_id=[];
    nn=length(varargin);
    for ii=1:2:nn
        if strcmp(varargin{ii},'shock_horizon')
            if nn>ii
                shock_horizon_id=[ii,ii+1];
            else
                shock_horizon_id=ii;
            end
        end
    end
    shock_horizon=varargin(shock_horizon_id);
    varargin(shock_horizon_id)=[];
    obj=set@rise_generic(obj,varargin{:});
    if ~isempty(shock_horizon_id)
        if numel(shock_horizon_id)==1
            % display the options
            error('somethings should be implemented in this case: contact junior.maih@gmail.com')
        else
            set_shock_horizon()
        end
    end
end

    function set_shock_horizon()
        value=shock_horizon{2};
        if isa(value,'double') && numel(value)==1
            value=round(abs(value));
            obj.exogenous.shock_horizon(~obj.exogenous.is_observed)=value;
        elseif isstruct(value)
            fields=fieldnames(value);
            shock_locs=locate_variables(fields,obj.exogenous.name);
            obs_status=obj.exogenous.is_observed(shock_locs);
            if any(obs_status)
                disp(fields(obs_status))
                error('the shocks above are deterministic and there horizon cannot be set')
            end
            for ifield=1:numel(fields)
                shock=fields{ifield};
                v=round(abs(value.(shock)));
                obj.exogenous.shock_horizon(shock_locs(ifield))=v;
            end
        else
            error('value must be a scalar or a structure')
        end
    end
end
