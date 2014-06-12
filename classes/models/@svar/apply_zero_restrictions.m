function [obj,linear_restrictions]=apply_zero_restrictions(obj)
if isempty(obj)
    obj=struct('restrict_lags',{{}},'estim_linear_restrictions',{{}});
    return
end
if ~strcmp(class(obj),'svar') %#ok<STISA>
    return
end
linear_restrictions=obj.options.restrict_lags;
param_template=obj.param_template;
if ~isempty(linear_restrictions)
    [nrest,ncols]=size(linear_restrictions);
    
    processed=false(nrest,1);
    
    % phase one apply the simple restrictions
    %----------------------------------------
    myconvert=@(z)int2str(locate_variables(z,obj.endogenous.name));
    operators={'+','-','*','/','^'};
    linear_restrictions(:,1)=svar.reformat_restriction(linear_restrictions(:,1),myconvert);
    for irest=1:nrest
        right=0;
        if ncols>1
            right=linear_restrictions{irest,2};
        end
        processed(irest)=apply_restriction(linear_restrictions{irest,1});
    end
    linear_restrictions=translate_svar_restrictions(linear_restrictions(~processed,:));
end
obj.estim_linear_restrictions=linear_restrictions;
obj.estim_param_template=param_template;

    function flag=apply_restriction(x)
        flag=false;
        % do not process if there is an operator
        %---------------------------------------
        for iop=1:length(operators)
            has_operator=any(x==operators{iop});
            if has_operator
                return
            end
        end
        [a_loc,eqtn_loc,var_loc]=svar.decompose_parameter(x,obj.param_template(1,:));
        param_template{2,a_loc}(eqtn_loc,var_loc)=right;
        flag=true;
    end

    function restrictions=translate_svar_restrictions(restrictions)
        % two things have to be done
        % 1- translate the remaining restrictions to form sparse matrices
        % of restrictions to be used during estimation
        names=obj.parameters.name;
        var_list=splanar.initialize(names,names);
        for irest_=1:size(restrictions,1)
            occur=regexp(restrictions{irest_,1},'\w+','match');
            args=ismember(names,occur);
            func=cell2mat(strcat(names(args),','));
            func=str2func(['@(',func(1:end-1),')',restrictions{irest_,1}]);
            args=var_list(args);
            zz=func(args{:});
            restrictions{irest_,1}=sparse(str2num(char(diff(zz,(1:obj.parameters.number))))); %#ok<ST2NM>
        end
    end

end
%         underscores=find(x=='_');
%         aa=x(1:underscores(1)-1);
%         a_loc= strcmp(aa,param_template(1,:));
%         eqtn_loc=str2double(x(underscores(1)+1:underscores(2)-1));
%         var_loc=str2double(x(underscores(2)+1:end));
% put in the form {lag,eqtn,var}
%         x=regexprep(x,'a(\w+)\((\d+),(\w+)\)','{$1,$2,${myconvert($3)}}');
%         x=regexprep(x,'a(\w+)\((\w+),(\w+)\)','{$1,${myconvert($2)},${myconvert($3)}}');