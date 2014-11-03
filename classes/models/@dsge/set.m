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
        elseif strcmp(varargin{ii},'solve_use_disc')
            use_disc_id=[ii,ii+1];
        end
    end
    solve_shock_horizon=varargin(shock_horizon_id);
    solve_use_disc=varargin(use_disc_id);
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
    if ~isempty(solve_use_disc)
        solve_use_disc=solve_use_disc{2};
        swap_routines();
    end
end

    function swap_routines()
        if isempty(obj.online_routines)
            obj.online_routines=obj.routines;
        end
        if solve_use_disc
            % do this systematically, irrespective of whether the routines
            % are the ones already written to disk or not
            write_routines_to_disk();
            obj.routines=obj.disc_routines;
        else
            obj.routines=obj.online_routines;
        end
        function write_routines_to_disk()
            MainFolder=obj.options.results_folder;
            obj.disc_routines=obj.routines;
            routines_names=fieldnames(obj.disc_routines);
            routines_names=setdiff(routines_names,{'symbolic','likelihood'});
            curr_dir=pwd();
            cd([MainFolder,filesep,'routines']);
            for irout=1:numel(routines_names)
                r_name=routines_names{irout};
                fname=r_name;
                rcode=utils.code.code2file(obj.disc_routines.(r_name),fname);
                if ~rcode
                    obj.disc_routines.(r_name)=str2func(['@',fname]);
                end
            end
            cd(curr_dir)
            % make the function readily available for use if necessary
            %---------------------------------------------------------
            rehash()
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
