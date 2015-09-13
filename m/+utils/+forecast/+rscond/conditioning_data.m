function [model,y0,ycond,econd,histdb,conddb]=conditioning_data(m,...
    y_names,e_names,end_history,do_soft)
if nargin<5
    do_soft=false;
end

% model solution
%----------------
model=model_solution(m);

% historical data
%-----------------
if ischar(y_names)
    y_names=cellstr(y_names);
end
if ischar(e_names)
    e_names=cellstr(e_names);
end
tight=2; % for the construction of bands
histdb=historical_data();

% initial conditions
%--------------------
y0=initial_condition();

% conditioning variables's positions
%-------------------------------------
[y_pos,e_pos]=condition_locations();

% observables and observables' data
%----------------------------------
cond_span=date2serial(end_history)+1:date2serial('2010Q3');
obsdb=histdb(y_names);
varnames=histdb.varnames;
conddb=obsdb;
if ~isempty(obsdb)
    conddb=obsdb(cond_span);
end
datay=build_cond_data(y_names);
datae=build_cond_data(e_names);
ycond=struct('data',datay,'pos',y_pos);
econd=struct('data',datae,'pos',e_pos);

    function datax=build_cond_data(names)
        nx=numel(names);
        nspan=numel(cond_span);
        datax=nan(nx,nspan,3);
        for ivar=1:nx
            vname=names{ivar};
            di=histdb(cond_span,vname);
            % always hard if not stated otherwise
            %------------------------------------
            datax(ivar,:,1)=di;
            datax(ivar,:,2)=di;
            datax(ivar,:,3)=di;
            % and possibly soft if possible
            %-------------------------------
            apply_soft(['lower_',vname],2)
            apply_soft(['upper_',vname],3)
            % should the user put -inf or inf themselves? most def but I
            % don't believe that if they put a number they will put
            % anything else than finite...
            % - Another interesting case if what happens if the user only
            % puts either the lower bound or just the upper bound.
            % - According to this specification, there is no conditionning
            % on a variable if the user does not give a central tendency.
            % - We have to revisit the role of nan central tendency vs nan
            % bounds...
        end
        function apply_soft(name,pos)
            low=any(strcmp(name,varnames));
            if low
                di_=histdb(cond_span,name);
                datax(ivar,:,pos)=di_;
            end
        end
    end

    function histdb=historical_data()
        m=filter(m);
        smooth_vals=m.filtering.Expected_smoothed_variables;
        smooth_shocks=m.filtering.Expected_smoothed_shocks;
        if do_soft
            for iname=1:numel(y_names)
                vname=y_names{iname};
                stdi=std(smooth_vals.(vname));
                smooth_vals.(['lower_',vname])=smooth_vals.(vname)-tight*stdi;
                smooth_vals.(['upper_',vname])=smooth_vals.(vname)+tight*stdi;
            end
            for iname=1:numel(e_names)
                vname=e_names{iname};
                % I might have many columns, pick one
                data_i=smooth_shocks.(vname)(:,1);
                stdi=std(data_i);
                smooth_shocks.(['lower_',vname])=data_i-tight*stdi;
                smooth_shocks.(['upper_',vname])=data_i+tight*stdi;
            end
        end
        histdb=utils.miscellaneous.mergestructures(smooth_vals,smooth_shocks);
        histdb=struct2pages(histdb);
        % For some reason I have many pages...
        histdb=histdb(:,:,1);
    end

    function [y,e]=condition_locations()
        y=model.inv_order_var(m.observables.state_id);
        ytmp=locate_variables(y_names,m.observables.name);
        y=y(ytmp);
        e=locate_variables(e_names,m.exogenous.name);
    end

    function y0=initial_condition()
        endo_names=m.endogenous.name(model.new_order);
        nvars=numel(endo_names);
        y0=nan(nvars,1);
        histvals=histdb(end_history);
        for ivar0=1:nvars
            vname=endo_names{ivar0};
            y0(ivar0)=double(histvals(vname));
        end
    end
end

    function model=model_solution(m)
        m=solve(m);
        [T,~,sstate,new_order,state_cols]=load_solution(m,'ov');
        iov=m.inv_order_var;
        Qfunc=prepare_transition_routine(m);
        % for the SVAR, there won't be a need to memoize
        Qfunc=memoize(Qfunc,iov);
        nshocks=sum(m.exogenous.number);
        model=struct('T',{T},'sstate',{sstate},'state_cols',state_cols,...
            'Q',m.solution.transition_matrices.Q,'Qfunc',Qfunc,...
            'k',max(m.exogenous.shock_horizon),...
            'new_order',new_order,'inv_order_var',iov,...
            'nshocks',nshocks);
    end