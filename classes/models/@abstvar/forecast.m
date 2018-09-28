function fkst=forecast(self,histdb,date_start,params,nsteps,...
    shock_uncertainty,Rfunc,conditions)

n=nargin;

set_defaults()

if ~isempty(histdb)

    self.estim_.data=histdb;

end

params=solve(self,params);

nreplic=numel(params);

[y,x,date_range_all]=collect_data(self,[]);

last_obs=find(date_range_all(1):date_range_all(2)==date_start-1);

ng=self.ng;

[y0,shocks]=format_input();

if ~isempty(conditions)

    yh=load_conditions();

end

for ig=1:ng

    if ig==1

        [fkst,info]=do_one_group(ig); %#ok<ASGLU>

        fkst=fkst(:,:,:,ones(1,ng));

    else

        fkst(:,:,:,ig)=do_one_group(ig);

    end

end

format_output()

    function [fkst,info]=do_one_group(ig)

        xig=[];

        if ~isempty(x)

            xig=x(:,:,ig);

        end

        if isempty(conditions)

            [fkst,info]=vartools.forecast(y0(:,:,ig),xig,params,shocks,Rfunc);

        else

            [fkst,info]=vartools.conditional_forecast(y0(:,:,ig),xig,params,Rfunc,yh(:,:,ig),shock_uncertainty);

        end

    end

    function yh=load_conditions()
        % the third dimension is for the panel case
        yh=nan(self.nvars,nsteps,ng);

        cnames=fieldnames(conditions);

        cond_span=date_start:date_start+nsteps-1;

        for iv=1:numel(cnames)

            vname=cnames{iv};

            crange=abstvar.reset_range(conditions.(vname));

            pos=strcmp(vname,self.endogenous);

            left=find(cond_span==crange(1));

            right=find(date_range_all(1):date_range_all(2)==crange(1));

            check_adequacy()

            nconstr=numel(crange(1):crange(2));
            % third dimension for the panel case
            yh(pos,left+(0:nconstr-1),:)=y(pos,right+(0:nconstr-1),:);

        end

        function check_adequacy()

            if crange(1)<cond_span(1)

                error('conditioning cannot start before the beginning of estimation')

            end

            if crange(2)>date_range_all(2)

                warning('conditioning cannot end after the available data: trimming from the end')

                crange(2)=date_range_all(2);

            end

            if crange(2)>cond_span(end)

                warning('conditioning cannot end after the forecast horizon: trimming from the end')

                crange(2)=cond_span(end);

            end

            if crange(2)<crange(1)

                error('end date of conditioning information cannot occur before start date')

            end

            if ~any(pos)

                error(['conditioning variable ',vname,' must be an endogenous variable'])

            end

        end

    end

    function format_output()

        fkst_out=struct();

        for g=1:ng

            tmp=struct();

            for iv=1:self.nvars

                if iv==1 && g==1

                    proto=ts(date_start-self.nlags,zeros(nsteps+self.nlags,nreplic));

                end

                dfkst=permute(fkst(iv,:,:,g),[2,3,1]);

                dfkst=[permute(y0(iv*ones(1,nreplic),:,g),[2,1,3]);dfkst]; %#ok<AGROW>

                tmp.(self.endogenous{iv})=set(proto,'data',dfkst);

            end

            if ng==1

                fkst_out=tmp;

            else

                fkst_out.(self.members{g})=tmp;

            end

        end

        fkst=fkst_out;

    end

    function [y0,shocks]=format_input()

        if ~isempty(x)

            xend=x(:,end,:);

            % there is a third dimension if we have panel...
            x=x(:,last_obs+1:end,:);

            ncols=size(x,2);

            % constant, if any, is already added by collect_data

            % expand the exogenous from the last observation if necessary
            if ncols<nsteps

                n_n=nsteps-ncols;

                xx=xend(:,ones(1,n_n),:);

                x=cat(2,x,xx);

                if ~all(x(:)==1)
                    % display warning only if it is not the constant
                    warning(['Not enough observations on exogenous over the ',...
                        'forecast horizon. Last observation replicated'])

                end

            end
            % there is a third dimension if we have panel...
            x=x(:,1:nsteps,:);

        end
        % there is a third dimension if we have panel...
        y0=y(:,1:last_obs,:);

        y0=y0(:,end-(self.nlags-1:-1:0),:);

        % there is a 4th dimension if we have panel...
        shocks=randn(self.nvars,nsteps,nreplic,ng);

        if ~shock_uncertainty

            shocks=zeros(size(shocks));

        end

    end

    function set_defaults()

        if n<8

            conditions=[];

            if n<7

                Rfunc=[];

                if n<6

                    shock_uncertainty=[];

                    if n<5

                        nsteps=[];

                        if n<4

                            params=[];

                            if n <3

                                date_start=[];

                                if n<2

                                    histdb=[];

                                end

                            end

                        end

                    end

                end

            end

        end

        % under unconditional forecast, there is not need for strong
        % identification
        if isempty(conditions),Rfunc=[]; end

        if isempty(Rfunc),Rfunc=identification(self,'choleski'); end

        if isempty(shock_uncertainty),shock_uncertainty=true; end

        if isempty(nsteps),nsteps=12; end

        if isempty(date_start),date_start=self.estim_.date_range(2)+1; end

        date_start=date2serial(date_start);

        if isempty(histdb),histdb=self.estim_.data; end

    end

end
