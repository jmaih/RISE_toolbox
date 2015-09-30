function obj=block_exogenous(endo_list,blk_exo_list,nlags)
% list of endogenous
% list of block exogenous
% number of lags
% this cannot be applied to markov switching. In that case, one
% has to explicitly list the coefficients that are zero. This
% so because different coefficients may be controled by
% different markov chains.
% To Do: extend to markov switching by exploiting the information on the
% parameters.
if ~isequal(sort(endo_list),endo_list)
    error('The list of endogenous must be sorted')
end
nv=numel(endo_list);
siz_obj=100;
incmnt=100;
obj=cell(siz_obj,2);
iter=0;
for ieqtn=1:nv
    foreign_equation=any(strcmp(endo_list{ieqtn},blk_exo_list));
    if ~foreign_equation
        continue
    end
    for ilag=1:nlags
        for ivar=1:nv
            endo_var=endo_list{ivar};
            home_variable=~any(strcmp(endo_var,blk_exo_list));
            if home_variable
                iter=iter+1;
                if iter==siz_obj
                    obj=[obj
                        cell(incmnt,2)]; %#ok<AGROW>
                end
                obj(iter,:)={coef(ieqtn,endo_var,ilag),0};
            end
        end
    end
end
% trim
obj=obj(1:iter,:);
end
