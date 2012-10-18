function model=load_parameters(model,the_mode_file)

% This loads the parameters and writes them to the parameter
% object. This allows the user to quickly load the parameters
% from a file, which may be the output of estimation, and get
% going with irfs, simulations, etc.
% This function eventually calls set_parameters, which contains
% alternative (online, i.e. without a file) ways of loading the
% parameters.
% There is another function that loads the parameters into
% the estimated parameters as start values. One nice feature
% is that it need not contain all the estimated parameters.

if isempty(model)
    model=struct();
    return
end
% read the mode file
fid=fopen(the_mode_file,'r');
markov_chains={model.markov_chains.name};
Regimes=model.Regimes;
params={'name','value','regime'};
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    tline(isspace(tline))=[];
    col=strfind(tline,':');
    par_val=str2double(tline(col+1:end));
    par_name=tline(1:col-1);
    par_chain='const';
    par_state=1;
    parenth=strfind(par_name,'(');
    if ~isempty(parenth)
        comma=strfind(par_name,',');
        item_1=par_name(parenth+1:comma-1);
        item_2=par_name(comma+1:end-1);
        par_name=par_name(1:parenth-1);
		par_chain=[]; par_state=[];
		if ismember(item_1,markov_chains)
			par_chain=item_1;
		else
			par_state=str2double(item_1);
		end 
		if ismember(item_2,markov_chains)
			par_chain=item_2;
		else
			par_state=str2double(item_2);
		end
		if isempty(par_state)||isempty(par_chain)
			error([mfilename,':: parameter information at line (',tline,') unrecognized'])
		end
    end
    % locate the regimes where par_name, controled by par_chain, is
    % in state par_state
    chain_loc=strcmp(par_chain,markov_chains);
    reg=find(Regimes(:,chain_loc)==par_state);
	if isempty(reg)
		error([mfilename,':: state ',int2str(par_state),' not in the range of states of markov chain ',par_chain])
	end
    for ireg=1:numel(reg)
        params=[params
            {par_name,par_val,reg(ireg)}];
    end
end
fclose(fid);

model=set_parameters(model,params);