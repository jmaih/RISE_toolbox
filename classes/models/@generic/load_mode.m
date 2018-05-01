function model=load_mode(model)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% reload the start values for the estimated parameters from a structure.

nobj=numel(model);
if nobj>1
    for iobj=1:nobj
        model(iobj)=load_mode(model(iobj));
    end
    return
end
the_mode_file=model.options.estim_start_vals;

if ~isempty(the_mode_file)
    if ~isstruct(the_mode_file)
        error('estim_start_vals should be a structure')
    end
	param_names={model.estimation.priors.name};
    % turn the official names into valid names. Since the mode file is a
    % structure, it has to be the case that it only includes valid names
    %---------------------------------------------------------------------
    param_names=parser.param_texname_to_param_name(param_names);
    re_started_names=fieldnames(the_mode_file);
    for irestart=1:numel(re_started_names)
        loc=find(strcmp(re_started_names{irestart},param_names));
        if isempty(loc)
            error(['parameter ',re_started_names{irestart},' is not an estimated parameter'])
        end
        model.estimation.priors(loc).start=the_mod_file.(re_started_names{irestart});
    end
end

% function model=load_mode(model)
% 
% the_mode_file=model.options.estim_mode_file;
% 
% if ~isempty(the_mode_file)&& exist(the_mode_file,'file') 
% 	% read the mode file
% 	fid=fopen(the_mode_file);
% 	param_names={model.estimated_parameters.name};
% 	while 1
% 	    tline = fgetl(fid);
% 	    if ~ischar(tline)
% 	        break
% 	    end
% 	    tline(isspace(tline))=[];
% 	    col=strfind(tline,':');
% 	    par_val=str2double(tline(col+1:end));
% 	    par_name=tline(1:col-1);
% 		par_loc=find(strcmp( par_name,param_names));
% 		if ~isempty(par_loc)
% 			model.estimated_parameters(par_loc).startval=par_val;
% 		end
% 	end
% 	fclose(fid);
% end