function model=load_mode(model)

the_mode_file=model.options.estim_mode_file;

if ~isempty(the_mode_file)&& exist(the_mode_file,'file') 
	% read the mode file
	fid=fopen(the_mode_file);
	param_names={model.estimated_parameters.name};
	while 1
	    tline = fgetl(fid);
	    if ~ischar(tline)
	        break
	    end
	    tline(isspace(tline))=[];
	    col=strfind(tline,':');
	    par_val=str2double(tline(col+1:end));
	    par_name=tline(1:col-1);
		par_loc=find(strcmp( par_name,param_names));
		if ~isempty(par_loc)
			model.estimated_parameters(par_loc).startval=par_val;
		end
	end
	fclose(fid);
end