function obj=translate_svar_restrictions(obj)
if isempty(obj)
    obj=struct();
    return
end
restr_name='lag_structure_zero_restrictions';
restrictions=obj.restrictions.(restr_name);
% translate the remaining restrictions to form sparse matrices
% of restrictions to be used during estimation
names=obj.parameters.name;
var_list=splanar.initialize(names,names);
for irest=1:size(restrictions,1)
    occur=regexp(restrictions{irest,1},'\w+','match');
    args=ismember(names,occur);
    func=cell2mat(strcat(names(args),','));
    func=str2func(['@(',func(1:end-1),')',restrictions{irest,1}]);
    args=var_list(args);
    zz=func(args{:});
    restrictions{irest,1}=sparse(str2num(char(diff(zz,(1:obj.parameters.number))))); %#ok<ST2NM>
end
obj.restrictions.(restr_name)=restrictions;
end