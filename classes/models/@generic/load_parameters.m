function model=load_parameters(model,the_mode_file)
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


% This loads the parameters. This allows the user to quickly
% load the parameters
% from a file, which may be the output of estimation, and get
% going with irfs, simulations, etc.


if isempty(model)
    
    model=cell(0,4);
    
    return
    
end
% read the mode file
fid=fopen(the_mode_file,'r');

params=struct();

while 1
    
    tline = fgetl(fid);
    
    if ~ischar(tline)
        
        break
        
    end
    
    tline(isspace(tline))=[];
    
    col=strfind(tline,':');
    
    par_val=str2double(tline(col+1:end));
    
    par_name=tline(1:col-1);
    
    % remove all blanks in the name
    par_name(isspace(par_name))=[];
    
    parenth=strfind(par_name,'(');
    
    if ~isempty(parenth)
        % what is represented in the model by alpha(coef,2)=10 then becomes
        % params.alpha_coef_2=10
        par_name=strrep(par_name,',','_');
        
        par_name=strrep(par_name,'(','_');
        
        par_name=strrep(par_name,')','');
        
    end
    
    params.(par_name)=par_val;
    
end

fclose(fid);

model=set(model,'parameters',params);

end