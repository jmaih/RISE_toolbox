function obj=set(obj,varargin)
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
% obj=set(obj,'solve_shock_horizon',struct('shock1',2,'shock3',4))
% obj=set(obj,'solve_shock_horizon',5)
%
% See also: 

if isempty(obj)
    obj=set@rise_generic(obj);
    % add other fields here if necessary
else
    % obj=set(obj,'solve_shock_horizon',struct('shock1',2,'shock3',4))
    % obj=set(obj,'solve_shock_horizon',5)
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
        swap_routines();
    end
end

    function swap_routines()
        routines_names=fieldnames(obj.online_routines);
        routines_names=setdiff(routines_names,{'symbolic','likelihood'});
        switch solve_function_mode
            case {'explicit','amateur'}
                obj.routines=obj.online_routines;
            case {'vectorized','professional'}
                vectorize_routines()
            case 'disc'
                write_routines_to_disk();
            otherwise
        end
        function vectorize_routines()
            for irout=1:numel(routines_names)
                r_name=routines_names{irout};
                obj.routines.(r_name)=utils.code.code2vector(obj.online_routines.(r_name));
            end
        end
        function write_routines_to_disk()
            MainFolder=obj.options.results_folder;
            routines_dir=[MainFolder,filesep,'routines'];
            for irout=1:numel(routines_names)
                r_name=routines_names{irout};
                fname=[r_name,'_',solve_function_mode,'__'];
                rcode=utils.code.code2file(obj.online_routines.(r_name),fname);
                if ~rcode
                    movefile([fname,'.m'],routines_dir,'f')
                    obj.routines.(r_name)=utils.code.func2fhandle([routines_dir,filesep,fname]);
                end
            end
        end
    end

    function set_shock_horizon()
        value=solve_shock_horizon{2};
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
