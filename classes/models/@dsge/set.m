function obj=set(obj,varargin)
% set - sets options for dsge|rise models
%
% Syntax
% -------
% ::
%
%   obj=set(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model object
%
% - **varargin** : valid input arguments coming in pairs. Notable fields to
%   that can be set include and are not restricted to:
%   - **solve_shock_horizon** [integer|struct|cell]
%       - for the integer case, all shocks are set to the same integer
%       - for the struct case, the input must be a structure with shock
%           names as fields. Only the shock names whose value is to change
%           have to be listed. In this case, different shocks can have
%           different horizons k. The default is k=0 i.e. agents don't
%           see into the future
%       - for the cell case, the cell should have two columns. The first
%           column includes the names of the shocks whose horizon is to
%           change. The second column includes the horizon for each shock
%           name on the left.
%   - **solve_function_mode** [{explicit/amateur}|vectorized/professional|disc]
%       - in the **amateur** or **explicit** mode the functions are kept in
%           cell arrays of anonymous functions and evaluated using for
%           loops
%       - in the **vectorized** or **professional** mode the functions are
%           compacted into one long and unreadable function.
%       - in the **disc** mode the functions are written to disc in a
%           subdirectory called routines.
%
% Outputs
% --------
%
% - **obj** [rise|dsge]: model object
%
% More About
% ------------
%
% Examples
% ---------
%
% obj=set(obj,'solve_shock_horizon',struct('shock1',2,'shock3',4))
% obj=set(obj,'solve_shock_horizon',5)
%
% See also: rise_generic.set

if isempty(obj)
    obj=set@rise_generic(obj);
    % add other fields here if necessary
    return
end

nobj=numel(obj);

shock_horizon_id=[];
use_disc_id=[];
nn=length(varargin);
for ii=1:2:nn
    if strcmp(varargin{ii},'solve_shock_horizon')
        if nn>ii
            shock_horizon_id=[ii,ii+1];
        else
            shock_horizon_id=ii;
        end
    elseif strcmp(varargin{ii},'solve_function_mode')
        use_disc_id=[ii,ii+1];
    end
end
solve_shock_horizon=varargin(shock_horizon_id);
solve_function_mode=varargin(use_disc_id);
varargin(shock_horizon_id)=[];
% do not remove the disc property since that option has to be visible
% unlike the shocks horizon. varargin(use_disc_id)=[];
obj=set@rise_generic(obj,varargin{:});
if ~isempty(shock_horizon_id)
    if numel(shock_horizon_id)==1
        % display the options
        error('somethings should be implemented in this case: contact junior.maih@gmail.com')
    else
        set_shock_horizon()
    end
end
if ~isempty(solve_function_mode)
    solve_function_mode=solve_function_mode{2};
    % do this one at a time
    %-----------------------
    for ii=1:nobj
        obj(ii)=swap_routines(obj(ii));
    end
end

    function this=swap_routines(this)
        routines_names=fieldnames(this.online_routines);
        routines_names=setdiff(routines_names,{'symbolic','likelihood'});
        switch solve_function_mode
            case {'explicit','amateur'}
                this.routines=this.online_routines;
            case {'vectorized','professional'}
                vectorize_routines()
            case 'disc'
                write_routines_to_disk();
            otherwise
        end
        function vectorize_routines()
            for irout=1:numel(routines_names)
                r_name=routines_names{irout};
                this.routines.(r_name)=utils.code.code2vector(this.online_routines.(r_name));
            end
        end
        function write_routines_to_disk()
            MainFolder=this.options.results_folder;
            routines_dir=[MainFolder,filesep,'routines'];
            for irout=1:numel(routines_names)
                r_name=routines_names{irout};
                fname=[r_name,'_',solve_function_mode,'__'];
                rcode=utils.code.code2file(this.online_routines.(r_name),fname);
                if ~rcode
                    movefile([fname,'.m'],routines_dir,'f')
                    this.routines.(r_name)=utils.code.func2fhandle([routines_dir,filesep,fname]);
                end
            end
        end
    end

    function set_shock_horizon()
        value=solve_shock_horizon{2};
        % turn into struct if cell
        %---------------------------
        process_cell()
        
        % this has to be done one at a time
        %-----------------------------------
        for iobj=1:nobj
            if isa(value,'double') && numel(value)==1
                value=round(abs(value));
                obj(iobj).exogenous.shock_horizon(~obj(iobj).exogenous.is_observed)=value;
            elseif isstruct(value)
                fields=fieldnames(value);
                shock_locs=locate_variables(fields,obj(iobj).exogenous.name);
                obs_status=obj(iobj).exogenous.is_observed(shock_locs);
                if any(obs_status)
                    disp(fields(obs_status))
                    error('the shocks above are deterministic and there horizon cannot be set')
                end
                for ifield=1:numel(fields)
                    shock=fields{ifield};
                    v=round(abs(value.(shock)));
                    obj(iobj).exogenous.shock_horizon(shock_locs(ifield))=v;
                end
            else
                error('value must be a scalar or a structure')
            end
        end
        
        function process_cell()
            if iscell(value)
                if size(value,2)~=2
                    error('when setting the shock horizon as a cell, the cell must have two columns')
                end
                ff=value(:,1);
                check=cellfun(@isvarname,ff,'uniformOutput',false);
                check=all([check{:}]);
                if ~check
                    error('when setting the shock horizon as a cell, the first column should have valid variable names')
                end
                value=cell2struct(value(:,2),ff,1);
            end
        end
    end
end
