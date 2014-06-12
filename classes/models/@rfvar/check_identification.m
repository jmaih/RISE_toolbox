function obj=check_identification(obj)
if isempty(obj)
    obj=struct();
    return
end
for iobj=1:numel(obj)
    obj(iobj)=set_structural_shocks(obj(iobj));
    obj(iobj)=translate_restrictions(obj(iobj));
    over_=false;
    under_=false;
    exact_=false;
    global_=false;
    if ~isempty(obj(iobj).nonlinear_restrictions)
        zero_restr=strcmp('lag_struct_Then_irf_zero_restr',{obj(iobj).nonlinear_restrictions.name});
        zero_restr=obj(iobj).nonlinear_restrictions(zero_restr);
        [Qordered,the_ranks]=rfvar.sort_Q(zero_restr.Q);
        the_ranks=the_ranks(:)';
        n=obj(iobj).endogenous.number;
        
        glob_flag=true;
        for jj = 1:n
            Qj = Qordered{1,jj};
            if glob_flag
                M=[Qj*zero_restr.f
                    eye(jj),zeros(jj,n-jj)];
                glob_flag=rank(M)==n;
            end
        end
        n_j=n-(1:n);
        over_=any(the_ranks - n_j > 0);
        under_=any(the_ranks - n_j < 0);
        exact_=all(the_ranks - n_j == 0);
        global_=glob_flag;
    end
    obj(iobj).identification=struct('over',over_,'exact',exact_,'under',under_,...
        'global',global_);
end
end
