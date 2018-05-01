function obj=set(obj,varargin)
% set - sets options for dsge|rise models
%
% ::
%
%
%   obj=set(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge]: model object
%
%    - **varargin** : valid input arguments coming in pairs. Notable fields to
%      that can be set include and are not restricted to:
%      - **solve_shock_horizon** [integer|struct|cell]
%          - for the integer case, all shocks are set to the same integer
%          - for the struct case, the input must be a structure with shock
%              names as fields. Only the shock names whose value is to change
%              have to be listed. In this case, different shocks can have
%              different horizons k. The default is k=0 i.e. agents don't
%              see into the future. The value of the fields is one of the
%              following:
%              - k
%              - {k,chainName_state}
%              - {chainName_state,k}
%          - for the cell case, the cell should have two columns. The first
%              column includes the names of the shocks whose horizon is to
%              change. The second column is one of the following:
%              - the horizon for the shock
%              - {horizon,chainName_state}
%              - {chainName_state,horizon}
%      - **solve_function_mode** [{explicit/amateur}|vectorized/professional|disc]
%          - in the **amateur** or **explicit** mode the functions are kept in
%              cell arrays of anonymous functions and evaluated using for
%              loops
%          - in the **vectorized** or **professional** mode the functions are
%              compacted into one long and unreadable function.
%          - in the **disc** mode the functions are written to disc in a
%              subdirectory called routines.
%
% Returns:
%    :
%
%    - **obj** [rise|dsge]: model object
%
% Note:
%
% Example:
%
%    obj=set(obj,'solve_shock_horizon',struct('shock1',2,'shock3',4))
%    obj=set(obj,'solve_shock_horizon',5)
%
%    See also: generic.set

if isempty(obj)
    
    obj=set@generic(obj);
    % add other fields here if necessary
    return
    
end

shock_horizon_id=[];

solve_log_approx_vars_id=[];

use_disc_id=[];

nn=length(varargin);

is_discard=false(1,nn);

for ii=1:2:nn
   
    this_option=varargin{ii};
    
    if strcmp(this_option,'solve_shock_horizon')
        
        if any([obj.is_dsge_var_model])
            
            error('Not possible to set "solve_shock_horizon" in a dsge-var model')
        
        end
        
        shock_horizon_id=[ii,ii+1];
        
        is_discard([ii,ii+1])=true;
    
    elseif strcmp(this_option,'solve_log_approx_vars') 
        
        if any([obj.is_dsge_var_model])
           
            error('Not possible to set "solve_log_approx_vars" in a dsge-var model')
        
        end
        
        solve_log_approx_vars_id=[ii,ii+1];
        
        is_discard([ii,ii+1])=true;
        
    elseif strcmp(this_option,'solve_function_mode')
        
        use_disc_id=[ii,ii+1];
        
    elseif any(strcmpi(this_option,{'priors','estim_priors'}))
        
        obj=setup_priors(obj,varargin{ii+1});
        
        is_discard([ii,ii+1])=true;
        
    end
    
end

solve_shock_horizon=varargin(shock_horizon_id);

solve_function_mode=varargin(use_disc_id);

solve_log_approx_vars=varargin(solve_log_approx_vars_id);

varargin(is_discard)=[];
% do not remove the disc property since that option has to be visible
% unlike the shocks horizon. varargin(use_disc_id)=[];

% do the following even if varargin isempty, it could well be the
% initialization phase
%--------------------------------------------------------------------------
obj=set@generic(obj,varargin{:});

nobj=numel(obj);

if ~isempty(shock_horizon_id)
    
    if numel(shock_horizon_id)==1
        % display the options
        
        error('something should be implemented in this case: contact junior.maih@gmail.com')
    
    else
        % turn into struct if cell
        %---------------------------
        value=process_cell(solve_shock_horizon{2});
        % do this one at a time
        %-----------------------
        for ii=1:nobj
            
            obj(ii)=set_shock_horizon(obj(ii));
        
        end
        
    end
    
end

if ~isempty(solve_log_approx_vars_id)
        % turn into struct if cell
        %---------------------------
        solve_log_approx_vars=solve_log_approx_vars{2};
        % do this one at a time
        %-----------------------
        for ii=1:nobj
            
            obj(ii)=set_log_approx_vars(obj(ii));
        
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
                
                if strcmp(r_name,'planner_osr_support')
                
                    continue
                
                end
                
                fname=[r_name,'_',solve_function_mode,'__'];
                
                rcode=utils.code.code2file(this.online_routines.(r_name),fname);
                
                if ~rcode
                
                    movefile([fname,'.m'],routines_dir,'f')
                    
                    this.routines.(r_name)=utils.code.func2fhandle([routines_dir,filesep,fname]);
                
                end
                
            end
            
        end
        
    end

    function this=set_log_approx_vars(this)
        
        this.endogenous.is_log_expanded(1,:)=false;
        
        pos=locate_variables(solve_log_approx_vars,this.endogenous.name);
        
        this.endogenous.is_log_expanded(pos)=true;
        
    end

    function this=set_shock_horizon(this)
        
        chain_names=this.markov_chains.chain_names;
        
        regimes=cell2mat(this.markov_chains.regimes(2:end,2:end));
        
        % this has to be done one at a time
        %-----------------------------------
        
        if isa(value,'double') && numel(value)==1
            
            value=round(abs(value));
            
            this.exogenous.shock_horizon(:,~this.exogenous.is_observed)=value;
        
        elseif isstruct(value)
            
            fields=fieldnames(value);
            
            shock_locs=locate_variables(fields,this.exogenous.name);
            
            obs_status=this.exogenous.is_observed(shock_locs);
            
            if any(obs_status)
                
                disp(fields(obs_status))
                
                error('the shocks above are deterministic and there horizon cannot be set')
            
            end
            
            for ifield=1:numel(fields)
                
                shock=fields{ifield};
                
                process_regime_horizon(value.(shock),shock_locs(ifield));
            
            end
            
        else
            
            error('value must be a scalar or a structure')
        
        end
        
        [~,this]=check_property(this,'solve_shock_horizon',value);
        
        function process_regime_horizon(vin,shockpos)
            
            if isa(vin,'double') && numel(vin)==1
                
                this.exogenous.shock_horizon(:,shockpos)=vin;
            
            elseif isa(vin,'cell') && size(vin,2)==2
                
                nrows=size(vin,1);
                
                for irow=1:nrows
                
                    if ischar(vin{irow,1})
                    
                        chain_state=vin{irow,1};
                    
                    elseif ischar(vin{irow,2})
                    
                        chain_state=vin{irow,2};
                    
                    end
                    
                    if isa(vin{irow,1},'double')
                    
                        v=vin{irow,1};
                    
                    elseif isa(vin{irow,2},'double')
                    
                        v=vin{irow,2};
                    
                    end
                    
                    chain_state1=regexprep(chain_state,'(\w+)\((\d+)\)','$1_$2');
                    
                    if ~ismember(chain_state1,this.markov_chains.state_names)
                    
                        error([chain_state,' is not recognized as a state'])
                    
                    end
                    
                    unders=find(chain_state1=='_');
                    
                    cn=chain_state1(1:unders-1);
                    
                    state=str2double(chain_state1(unders+1:end));
                    
                    chainloc=strcmp(cn,chain_names);
                    
                    state_pos= regimes(:,chainloc)==state;
                    
                    this.exogenous.shock_horizon(state_pos,shockpos)=v;
                
                end
                
            else
                
                disp(vin)
                
                error('wrong format for the shock horizon')
            
            end
            
        end
        
    end

end

function value=process_cell(value)

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
