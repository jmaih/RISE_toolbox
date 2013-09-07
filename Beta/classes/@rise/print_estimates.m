function print_estimates(model,type,file2save2)
if nargin<3
    file2save2=[];
    if nargin<2
        type=[];
    end
end
precision=9;
if isempty(type)
    type='mode';
end
if ~ismember(type,{'mode','mean','post_sim_mode','prior_mean','startval'})
    error([mfilename,':: unrecognized type ',type])
end

param_names=char({model.estimated_parameters.name}');
param_vals=num2str(vertcat(model.estimated_parameters.(type)),precision);
params=[param_names,repmat(' : ',size(param_vals,1),1),param_vals];
% params=strcat(param_names,' : ',param_vals); This will not preserve
% trailing blanks
if ~isempty(file2save2)
    ext_loc=strfind(file2save2,'.');
    extension='out';
    if ~isempty(ext_loc)
        extension=file2save2(ext_loc+1:end);
        file2save2=file2save2(1:ext_loc-1);
    end
    if length(file2save2)>5
        if ~strcmp(file2save2(end-(4:-1:0)),'_mode')
            file2save2=[file2save2,'_mode'];
        end
    end
    fid=fopen([file2save2,'.',extension],'w');
else
    fid=1;
end

for ii=1:size(params,1)
    fprintf(fid,'%s \n',params(ii,:));
end
if ~isempty(file2save2)
    fclose(fid);
end
