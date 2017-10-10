function [Rfunc,ident]=identification(nvars,nlags,nx,restrictions,agnostic,max_trials,debug)

n=nargin;

set_defaults()

quick_exit=false;

restrictions=restrictions2restrictions(restrictions);

number_of_zero_restrictions=0;

[is,loc]=find_positions();

number_of_sign_restrictions=sum(is.sign);

rzlag=[]; rslag=[]; rllag=[];

collect_lags_for_restrictions()

rlag=[unique([rslag,rzlag]),unique(rllag)];

reset_lags()

nrows_phibar=nvars*numel(rlag); % nrows_phibar=size(PHIBAR,1)+nvars*numel(PHIBARzr);

[ZZ,itags_zero,ident]=build_z_matrices();

finite=@(x)x(isfinite(x) & imag(x)==0);

order=max(finite(rlag));

Rfunc=@identification_engine;

    function [RQP,retcode]=identification_engine(p)

		B=p.B;

		SIG=p.S;
        
        PHI=vartools.moving_average(B,nx,order);
        
        [PHIBAR,PHIBARzr]=add_restrictions_impact();
        
        trials=0;
        
        check=false;
        
        R=chol(SIG,'lower');
        
        while ~check && trials<max_trials
            
            trials=trials+1;
            
            RQ=generate_rotation();
            
            P=var_rotation(RQ);
            
            % by construction, zero restrictions are strictly enforced
            %---------------------------------------------------------
            RQP = RQ*P;
            
            % check sign restrictions
            %------------------------
            check=process_sign_restrictions(RQP);
            
        end
        
        retcode=1-check;
        
        if debug
            
            disp([trials,check])
            
            disp('R*R-SIG')
            
            max(max(abs(RQP*RQP.'-SIG)))
            
        end
        
        function PHIB=aggregate(R)
            
            PHIB=[];
            
            if ~isempty(PHIBAR)
                
                PHIB=PHIBAR*R;
                
            end
            
            if ~isempty(PHIBARzr)
                
                Ri=R\eye(nvars);
                
                tmp=PHIBARzr;
                
                for ii=1:numel(PHIBARzr)
                    
                    tmp{ii}=transpose(Ri*tmp{ii});
                    
                end
                
                PHIB=[PHIB;cell2mat(tmp)];
                
            end
            
        end
        
        function [PHIBAR,PHIBARzr]=add_restrictions_impact()
            
            PHIBAR=add_shock_restrictions();
            
            PHIBARzr=add_lag_structure_zero_restrictions();
            
            function PHIBAR=add_shock_restrictions()
                
                PHIBAR=cell(order+1,1);
                
                iterz=0;
                
                for ilag=1:order+1 % starting from 0
                    
                    if ismember(ilag-1,rlag)
                        
                        iterz=iterz+1;
                        
                        PHIBAR{iterz}=PHI{ilag};
                        
                    end
                    
                end
                
                PHIBAR=cell2mat(PHIBAR(1:iterz));
                
                if ismember(inf,rlag)
                    
                    Xinf=vartools.long_run_impact(B,eye(nvars),nx);
                    
                    PHIBAR=[PHIBAR;Xinf];
                    
                end
                
            end
            
            function PHIBARzr=add_lag_structure_zero_restrictions()
                
                if isempty(rllag)
                    
                    PHIBARzr=[];
                    
                    return
                    
                end
                
                % restrictions can be imposed on Cx, A0, A1,...,Ap but imposing
                % restrictions on Cx does not contribute to identification and as a
                % result we do no implement that just yet. Therefore we discard the
                % deterministic columns
                %
                % add the contemporaneous batch in Ci
                IB=[eye(nvars),B(:,nx+1:end)];
                
                if real(rllag(end))>nlags
                    
                    error(['Not possible to impose lag structure ',...
                        'restrictions on lags beyond ',int2str(nlags)])
                    
                end
                
                PHIBARzr=cell(nlags+1,1);
                
                iter=0;
                
                for ilag=1:nlags+1 % starting from 0
                    
                    if ismember(ilag-1+1i,rllag)
                        
                        iter=iter+1;
                        
                        PHIBARzr{iter}=IB(:,(ilag-1)*nvars+1:ilag*nvars);
                        
                    end
                    
                end
                
                PHIBARzr=PHIBARzr(1:iter);
                
            end
            
        end
        
        function check=process_sign_restrictions(RQP)
            
            check=false(1,number_of_sign_restrictions);
            
            if number_of_sign_restrictions
                
                if ~(ident.isunder||ident.isexact)
                    
                    error(['The model must be under or exactly identified '...
                        'to be able to do sign restrictions'])
                    
                end
                
                PHITILDE=aggregate(RQP);
                
                for irest=1:number_of_sign_restrictions
                    
                    signrestr=restrictions(loc.sign(irest));
                    
                    S=zeros(1,nrows_phibar);
                    
                    S((signrestr.lag_image-1)*nvars+signrestr.vbl)=1;
                    
                    issign=numel(signrestr.rhs)==1;
                    
                    SPHI=S*PHITILDE(:,signrestr(1).shk);
                    
                    if issign
                        
                        check(irest)=sign(SPHI)==sign(signrestr.rhs);
                        
                    else
                        
                        check(irest)=SPHI>=signrestr.rhs(1) &&...
                            SPHI<=signrestr.rhs(2);
                        
                    end
                    
                    if ~check(irest)
                        
                        break
                        
                    end
                    
                    if debug
                        
                        %                     disp([trials,irest,SPHI])
                        
                    end
                    
                end
                
            end
            
            check=all(check);
            
        end
        
        function P=var_rotation(RQ)
            
            if number_of_zero_restrictions
                
                PHIBAR_RQ=aggregate(RQ);
                
                P=var_zero_rotation(PHIBAR_RQ,ZZ,nvars,itags_zero,agnostic,quick_exit);
                
            else
                
                P=eye(nvars);
                
            end
            
        end
        
        function RQ=generate_rotation()
            
            % random matrix
            X=randn(nvars); % <--- X=rand(nvars);
            
            [QQQ,RR]=qr(X);
            
            for ii = 1:nvars
                
                if RR(ii,ii)<0
                    
                    QQQ(:,ii) = -QQQ(:,ii);
                    
                end
                
            end
            
            RQ=R*QQQ;
            
        end
        
    end

    function [Z,itags_zero,ident]=build_z_matrices()
        
        Z=struct();
        
        for ishock=1:nvars
            
            iter=0;
            
            Zi=zeros(nvars,nrows_phibar);
            
            % pure zero restrictions
            %-----------------------
            if any(is.zero)
                
                if ishock==1
                    
                    rzshocks=[restrictions(is.zero).shk];
                    
                end
                
                assign_shock(restrictions(is.zero),rzshocks)
                
            end
            
            % lag structure restrictions
            %---------------------------
            if any(is.lagstr)
                
                if ishock==1
                    
                    rlshocks=[restrictions(is.lagstr).shk];
                    
                end
                
                assign_shock(restrictions(is.lagstr),rlshocks)
                
            end
            
            Zi=Zi(1:iter,:);
            
            Z(ishock).rank=rank(Zi);
            
            Z(ishock).Z=sparse(Zi);
            
        end
        
        [Z,itags_zero]=sort_Z(Z);
        
        ident=check_identification([Z.rank],nvars);
        
        function assign_shock(rs,rsshocks)
            
            rsi=find(rsshocks==ishock);
            
            if ~isempty(rsi)
                
                for jj=1:numel(rsi)
                    
                    iter=iter+1;
                    
                    rsij=rs(rsi(jj));
                    
                    Zi(iter,(rsij.lag_image-1)*nvars+rsij.vbl)=1;
                    
                end
                
            end
            
        end
        
    end

    function collect_lags_for_restrictions()
        
        % pure zero restrictions
        %------------------------
        rzlag=[restrictions(is.zero).lag];
        
        number_of_zero_restrictions=number_of_zero_restrictions+...
            sum(is.zero);
        
        % pure sign restrictions
        %------------------------
        rslag=[restrictions(is.sign).lag];
        
        % adding lag structure in the form of complex numbers
        %-----------------------------------------------------
        % add a complex to every lag restriction...
        %------------------------------------------
        rllag=[restrictions(is.lagstr).lag]+1i;
        tmp=num2cell(rllag);
        % this is slick!!!
        %-----------------
        [restrictions(is.lagstr).lag]=(tmp{:});

        number_of_zero_restrictions=number_of_zero_restrictions+...
            sum(is.lagstr);
        
    end

    function reset_lags()
        
        % force a new field
        restrictions(1).lag_image=nan;
        
        restrictions(is.sign)=reset_lag_engine(restrictions(is.sign),rlag);
        
        restrictions(is.zero)=reset_lag_engine(restrictions(is.zero),rlag);
        
        restrictions(is.lagstr)=reset_lag_engine(restrictions(is.lagstr),rlag);
        
    end

    function set_defaults()
        
        if n<7
            
            debug=[];
            
            if n<6
                
                max_trials=[];
                
                if n<5
                    
                    agnostic=[];
                    
                    if n<4
                        
                        restrictions=[];
                        
                    end
                    
                end
                
            end
            
        end
        
        if isempty(debug),debug=false; end
        
        if isempty(agnostic),agnostic=true; end
        
        if isempty(max_trials),max_trials=10000; end
        
    end

    function [is,loc]=find_positions()
        
        types={restrictions.type};
        
        is.sign=strcmp(types,'s')|strcmp(types,'m'); % sign and magnitude
        
        is.zero=strcmp(types,'z'); % zero restrictions
        
        is.lagstr=strcmp(types,'ls'); % zero restrictions
        
        loc.zero=find(is.zero);
        
        loc.sign=find(is.sign);
        
        loc.lagstr=find(is.lagstr);
        
    end

end

function r=restrictions2restrictions(r)

if isempty(r)
    
    return
    
end

if ~isstruct(r)
    
    error('restrictions must be a structure')
    
end

r=decentralize(r);

    function r1=decentralize(r0)
        
        n=numel(r0);
        
        r1=r0(1:0); % initialize in a way that preserves the field names
        
        iter=0;
        
        for kk=1:n
            
            lags=r0(kk).lag;
            
            tmp=r0(kk);
            
            for ll=1:numel(lags)
                
                tmp.lag=lags(ll);
                
                iter=iter+1;
                
                r1(iter)=tmp;
                
            end
            
        end
        
    end


end

function P=var_zero_rotation(fa0ap,Z,nvars,itags_zero,agnostic,quick_exit)

P=zeros(nvars);

for jj=1:nvars
    
    Zj=Z(jj).Z;
    
    if quick_exit && isempty(Zj)
        
        warning('Leaving columns of zeros in the rotation is not a good idea')
        
        break
        
    end
    
    Mj=[Zj*fa0ap
        P(:,1:jj-1)'];
    
    P(:,jj)=find_vector();
        
end

% re-order because the order of ZZ was re-ordered
P=P(:,itags_zero);

if any(isnan(P(:)))
    
    error('nans in the zero rotation matrix. Maybe something wrong with the restrictions.')
    
end

% fa0ap*P

    function v=find_vector()
        
        if agnostic
            
            Nj=null(Mj);
            
            xj=randn(nvars,1);
            
            Njxj=Nj.'*xj;
            
            NjNjxj=Nj*Njxj/norm(Njxj);
            
            v=NjNjxj;
            
        else
            
            [q,~]=qr(Mj');
            
            v=q(:,end);
            
        end
        
    end

end

function ident=check_identification(the_ranks,nvars)

batch=the_ranks(:).'-(nvars-1:-1:0);

ident=struct();

ident.isexact=all(batch==0);

ident.isover=any(batch>0);

ident.isunder=any(batch<0);

end

function rr=reset_lag_engine(rr,batch)

if isempty(rr)
    
    return
    
end

rrlags=[rr.lag];

if ~all(ismember(rrlags,batch))
    
    error('lag not found in the set of lags')
    
end
    
rrlags_image=rrlags;

for iii=1:length(batch)
    
    ping=rrlags==batch(iii);
    
    rrlags_image(ping)=iii;
    
end

rrlags_image=num2cell(rrlags_image);

% this is slick!!!
%-----------------
[rr.lag_image]=(rrlags_image{:});

end

function [Z,itags]=sort_Z(Z)

the_ranks=[Z.rank];

[~,tags]=sort(the_ranks(:),1,'descend');

Z=Z(tags);

itags(tags)=1:numel(tags);

end

%{


function r=restrictions2restrictions(r)

if isempty(r)
    
    return
    
end

if ~isstruct(r)
    
    error('restrictions must be a structure')
    
end

ff=fieldnames(r);

default_fields={'zeros','sign','lag_structure'};

zero_fields={'v','shock','lag'}; % order not to be modified

sign_fields=[zero_fields,'impact']; % order not to be modified

lag_struct_fields={'v','eqtn','lag'}; % order not to be modified

if ~all(ismember(ff,default_fields))
    
    error('fields of restrictions must be members of ''zeros'',''sign'',''lag_structure''')
    
end

for ii=1:numel(default_fields)
    
    if ~isfield(r,default_fields{ii})
        
        r.(default_fields{ii})=[];
        
    end
    
    if isempty(r.(default_fields{ii}))
        
        continue
        
    end
    
    if iscell(r.(default_fields{ii}))
        
        do_cell2struct(default_fields{ii})
        
    elseif isstruct(r.(default_fields{ii}))
        
        chkstruct(default_fields{ii})
        
    else
        
        error([default_fields{ii},' field expected to be a cell or a struct'])
        
    end
    
    % for sign restrictions replace + and - with 1 and -1
    %----------------------------------------------------
    if strcmp(default_fields{ii},'sign')
        
        for jj=1:numel(r.sign)
            
            impact=r.sign(jj).impact;
            
            if ischar(impact)
                
                r.sign(jj).impact=if_elseif(strcmp(impact,'+'),1,...
                    strcmp(impact,'-'),-1);
            else
                % magnitude restrictions
                r.sign(jj).impact=sort(impact);
                
            end
            
        end
        
    end
    
    % decentralize in case there are many lags in one restriction
    %------------------------------------------------------------
    r.(default_fields{ii})=decentralize(r.(default_fields{ii}));
    
end

    function r1=decentralize(r0)
        
        n=numel(r0);
        
        r1=r0(1:0); % initialize in a way that preserves the field names
        
        iter=0;
        
        for kk=1:n
            
            lags=r0(kk).lag;
            
            tmp=r0(kk);
            
            for ll=1:numel(lags)
                
                tmp.lag=lags(ll);
                
                iter=iter+1;
                
                r1(iter)=tmp;
                
            end
            
        end
        
    end

    function do_cell2struct(sname)
        
        switch sname
            
            case  'zeros'
                
                target_fields=zero_fields;
                
            case 'sign'
                
                target_fields=sign_fields;
                                
            case 'lag_structure'
                
                target_fields=lag_struct_fields;
                
        end
        
        s=r.(sname);
        
        [~,ncols]=size(s);
        
        if ncols~=numel(target_fields)
            
            error(['wrong number of columns for cell ',sname])
            
        end
        
        r.(sname)=cell2struct(s,target_fields,2);
        
    end

    function chkstruct(sname)
        
        s=r.(sname);
        
        ffs=sort(fieldnames(s));
        
        switch sname
            
            case 'zeros'
                
                if ~all(strcmp(ffs(:).',sort(zero_fields)))
                    
                    error(['valid fields for "',sname,'" restrictions are ''v'',''shock'',''lag'''])
                    
                end

            case 'sign'
                
                if ~all(strcmp(ffs(:).',sort(sign_fields)))
                    
                    error(['valid fields for "',sname,'" restrictions are ''v'',''shock'',''lag'',''impact'''])
                    
                end
                
            case 'lag_structure'
                
                if ~all(strcmp(ffs(:).',sort(lag_struct_fields)))
                    
                    error(['valid fields for "',sname,'" restrictions are ''v'',''eqtn'',''lag'''])
                    
                end
                
        end
        
    end

end

%}