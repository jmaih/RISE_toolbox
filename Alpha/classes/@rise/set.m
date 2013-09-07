function obj=set(obj,property,value)
if isempty(obj)
    obj=struct();
    return
end
nobj=numel(obj);
if nobj>1
    for iobj=1:nobj
        obj(iobj)=set(obj(iobj),property,value);
    end
    return
end
if ismember(property,{'state_tex_names','chain_tex_names','regime_tex_names'})
    if ischar(value)
        value=cellstr(value);
    end
    name=obj.markov_chains.(strrep(property,'tex_',''));
    if numel(name)~=numel(value)
        error(['number of elements in ',property,' should be the same as the number of elements in ',name])
    end
    obj.markov_chains.(property)=value(:)';
elseif strcmpi(property,'parameters')
    obj=setup_calibration(obj,value);
%     fields=fieldnames(value);
%     parameter_values=obj.parameter_values;
%     par_list=get(obj,'par_list');
%     chain_list=obj.markov_chains.chain_names;
%     governing_chain=obj.parameters.governing_chain;
%     regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));
%     for ifield=1:numel(fields)
%         pname0=fields{ifield};
%         pname=pname0;
%         loc=find(strcmp(pname,par_list));
%         chain='const';
%         state=1;
%         if isempty(loc)
%             underscores=strfind(pname,'_');
%             if isempty(underscores)||numel(underscores)<2
%                 error(['parameter ''',pname,''' could not be located'])
%             end
%             state=str2double(pname(underscores(end)+1:end));
%             chain=pname(underscores(end-1)+1:underscores(end)-1);
%             pname=pname(1:underscores(end-1)-1);
%             loc=find(strcmp(pname,par_list));
%         end
%         chain_position=governing_chain(loc);
%         if ~strcmp(chain,chain_list{chain_position})
%             error(['parameter ''',pname,''' is not controlled by markov chain ''',chain,''])
%         end
%         destination=regimes(:,chain_position)==state;
%         parameter_values(loc,destination)=value.(pname0);
%     end
%     obj.parameter_values=parameter_values;
else
    error(['',property,''' is not a settable property of RISE'])
end
end