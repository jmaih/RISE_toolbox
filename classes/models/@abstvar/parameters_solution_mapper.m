function [out,markov_chains]=parameters_solution_mapper(p,mc_chains,...
    endo_names,exo_names,panel_members,nx)

% N.B: nx =numel(exo_names)+constant(s) : in a panel, there could be
% multiple "constants" !!!
%
% for panel, push the endo_names followed by the country names to be
% consistent with the way the data are built...
%
% then if there is a restriction that contains an endogenous variable name
% it should be multiplicated
%
% what will not be allowed is for the restrictions across different panel
% units

np=numel(panel_members);

old_endo_names=endo_names;

old_exo_names=exo_names;

all_names=recell2matize([old_endo_names(:).',old_exo_names(:).']);

recreate_names()

% a, b
% a1, b2
% a(1,:), a(1)
% a(1,3)
% a(1,gdp)
% a(:,gdp)
% a2(1,gdp)
% a2(:,gdp)
% c
% c(1)
% c(1,ffr)
% c(:,ffr)
% s
% s(1),s(1,:)
% s(1,2)
% check whether there is a constant

out=struct();

pnames=p.pnames;

% cc=regexp(pnames,'c_\d+_0','match'); cc=[cc{:}]; constant=~isempty(cc);

% check whether it is a svar
% cc=regexp(pnames,'s_\d+_\d+','match'); cc=[cc{:}];
% 
% issvar=isempty(cc);

issvar=any(strcmp(pnames,'a0_1_1'));

rep_ab=@replace_ab_params; %#ok<NASGU>

rep_c=@replace_c_params; %#ok<NASGU>

rep_s=@replace_s_params; %#ok<NASGU>

edn=recell2matize(endo_names);

exn=recell2matize(exo_names);

if ~isempty(exn),exn=[exn,'|'];end

okura='\d+(\:\d+)?|\[[\d+,]+\]|\:';

% a, a8, a8(:),a8(3),a8(3,:),a8(3,C),a8([3,1,5]),a8(3:5,C)....
expressParamsAB=['\<(a|b)(\d+)?\>(?:\()?(',okura,')?,?(',edn,...
    '|\d+|\:)?(?:\))?'];

expressParamsC=['\<(c)\>(?:\()?(',okura,')?,?(',exn,'\d+|\:)?(?:\))?'];

if issvar
    
    expressParamsS=['\<(s)\>(?:\()?(',okura,')?(?:\))?'];
    
else
    
    expressParamsS=['\<(s)\>(?:\()?(',okura,')?,?(',okura,')?(?:\))?'];
    
end

nc=numel(mc_chains);

pnames=pnames(:).';

out.nparams=numel(pnames);

isTaken=false(1,out.nparams);

markov_chains=struct.empty;

pockets=group_parameters();

% After grouping sort the chains
cnames_={markov_chains.name};

pos=locate_variables(sort(cnames_),cnames_);

markov_chains=markov_chains(pos);

[out.regimes,out.journal,out.nregimes]=create_regimes();

% if any pocket is empty scream
bad=cellfun(@isempty,pockets(2:end));

if any(bad)
    
    badnames={mc_chains(bad).name};
    
    if ~isempty(badnames)
        
        disp(badnames)
        
        error('the chain above do not control any parameter')
        
    end

end

out.solution=do_solution(p,markov_chains,...
    mc_chains,...
    issvar,numel(endo_names),nx);

    function recreate_names()
        
        if np==0
            
            return
            
        end
        
        endo_names=do_one_type(endo_names);
        
        exo_names=do_one_type(exo_names);
        
        function vout=do_one_type(vin)
            
            nv=numel(vin);
            
            nvp=np*nv;
            
            vout=cell(1,nvp);
            
            iter=0;
            
            for iv=1:nv
                
                vv=vin{iv};
                
                for ip=1:np
                    
                    iter=iter+1;
                    
                    vout{iter}=[vv,'_',panel_members{ip}];
                
                end
                
            end
        
        end
        
    end

    function [regimes,journal,nregimes]=create_regimes()
                
        cnames={markov_chains.name};
        
        v=[markov_chains.number_of_states];
        
        nregimes=prod(v);
        
        regimes=cell(nregimes+1,nc+2);
        
        regimes(1,2:end)=cnames;
        
        for ireg=1:nregimes
            
            regimes{ireg+1,1}=sprintf('regime%0.0f',ireg);
            
        end
        
        [R,journal]=utils.gridfuncs.chain_grid(v);
        
        regimes(2:end,2:end)=num2cell(R);
        
    end

    function tmp=extract_parameters(cp)
        
        if ischar(cp)
            
            cp=cellstr(cp);
            
        end
        
        % transform
        cp=regexprep(cp,expressParamsAB,'${rep_ab($1,$2,$3,$4)}');
        
        cp=regexprep(cp,expressParamsC,'${rep_c($1,$2,$3)}');
        
        if issvar
            
            cp=regexprep(cp,expressParamsS,'${rep_s($1,$2)}');
            
        else
            
            cp=regexprep(cp,expressParamsS,'${rep_s($1,$2,$3)}');
            
        end
        
        tmp=cell(1,0);
        
        cp=expand_code(cp);
        
        for ii=1:numel(cp)
            
            extr=regexp(pnames,['\<',cp{ii},'\>'],'match');
            
            extr=[extr{:}];
            
            if isempty(extr)
                
                error('parameters not found or expression could not be parsed successfully')
                
            end
            
            locs=locate_variables(extr,pnames);
            
            if any(isTaken(locs))
                
                error('parameters declared multiple times or controlled by different mc chains')
                
            end
            
            isTaken(locs)=true;
            
            tmp=[tmp,extr(:).']; %#ok<AGROW>
            
        end
        
    end

    function out=replace_s_params(s,eqtn,col)
        
        if nargin<3,col=[]; end
        
        if issvar
            
            if isempty(eqtn)
                
                if ~isempty(col)
                    
                    error('s: eqtn empty and col not empty')
                    
                end
                
            else
                
                if isempty(col)
                    
                    col=eqtn;
                    
                elseif ~strcmp(eqtn,col)
                    
                    error('s: eqtn must be the same as col')
                    
                end
                
            end
            
        else
            
            if ~((isempty(col) && isempty(eqtn))||...
                    (~isempty(col) && ~isempty(eqtn)))
                
                error(['for reduced-form VARs, one should have both indices ',...
                    'either empty or non-empty'])
                
            end
            
            
        end
        
        dcol=~isempty(col) && all(isstrprop(col,'digit'));
        
        deqtn=~isempty(eqtn) && all(isstrprop(eqtn,'digit'));
        
        success=try_swap(); % change s_i_j into s_j_i if i<j
        
        eqtnx=update_equation(eqtn);
        
        colx=update_equation(col);
        
        out=[s,'_',eqtnx,'_',colx];
        
        if success||all(strcmp({eqtnx,colx},'\d+'))
            
            % do nothing
            
        elseif (dcol && strcmp(eqtnx,'\d+'))||...
                (deqtn && strcmp(colx,'\d+'))
            
            out=['(',out,'|',s,'_',colx,'_',eqtnx,')'];
            
        else
            
            error('unable to successfully parse expression')
            
        end
        
        
        function success=try_swap()
            % covariances are stored in vech form and so we may miss s_i_j
            % if j<i
            
            success=false;
            
            if dcol && deqtn
                
                tmpcol=str2double(col);
                
                tmpeqtn=str2double(eqtn);
                
                if tmpcol>=tmpeqtn
                    
                    col=int2str(tmpeqtn);
                    
                    eqtn=int2str(tmpcol);
                    
                    success=true;
                    
                end
                
            end
            
        end
        
    end

    function out=replace_c_params(c,eqtn,v)
        
        vx=update_variable(v,exo_names);
        
        % constant is always referred to as 0
                
        eqtnx=update_equation(eqtn);
        
        out=[c,'_',eqtnx,'_',vx];
        
    end

    function out=replace_ab_params(ab,lag,eqtn,v)
        
        lagx='\d+';
        
        update_lag()
        
        vx=update_variable(v,endo_names);
        
        eqtnx=update_equation(eqtn);
        
        out=[ab,lagx,'_',eqtnx,'_',vx];
        
        function update_lag()
            
            if ~isempty(lag)
                
                if strcmp(lag,':')
                    
                    % don't do anything
                    
                else
                    
                    lagx=lag;
                    
                end
                
            end
            
        end
        
    end

    function pockets=group_parameters()
        
        pockets=cell(1,nc+1);
        
        for ii=1:nc
            
            cp=mc_chains(ii).controlled_parameters;
            
            cp=big_extraction(cp);
            
            pockets{ii+1}=cp(:);
            
            set_chain()
            
        end
        
        % non-switching parameters
        pockets{1}=pnames(~isTaken);
        
        const_mc=parser.initialize_markov_chain('const',1,'is_endogenous',false,...
                'param_list',pockets{1}(:).','is_switching',false,...
                'duration',inf);
        
        markov_chains=[
            const_mc,...
                markov_chains
                ];
            
        function xout=big_extraction(xin)
                
                nxx=numel(xin);
                
                xout=cell(1,nxx);
                
                for jj=1:nxx
                    
                    xout{jj}=xpand(xin{jj});
                    
                end
                
                xout=[xout{:}];
                
                xout=extract_parameters(xout);
                
            function x=xpand(x)
                
                if ischar(x)
                    
                    x=cellstr(x);
                    
                end
                
                if np>1
                    
                    x=xpand_variables(x);
                    
                    x=decellize(x);
                    
                    x=xpand_digits(x);
                    
                    x=decellize(x);
                    
                end
                
                function x=xpand_digits(x)
                    
                    nxxx=numel(x);
                    
                    outx=cell(1,nxxx);
                    
                    for jjj=1:nxxx
                        
                        [ss,ee]=regexp(x{jjj},'\(\d+','start','end');
                        
                        if isempty(ss)
                            
                            outx{jjj}=x(jjj);
                            
                        else
                            
                            plg=x{jjj}(1:ss);
                            
                            mdl=str2double(x{jjj}(ss+1:ee));
                            
                            mdl=(mdl-1)*np+1:mdl*np;
                            
                            eplg=x{jjj}(ee+1:end);
                            
                            tmp=cell(1,np);
                            
                            for ip=1:np
                                
                                tmp{ip}=sprintf('%s%0.0g%s',plg,mdl(ip),eplg);
                                
                            end
                            
                            outx{jjj}=tmp;
                            
                        end
                    
                    end
                    
                    x=[outx{:}];
                    
                end
                
                function x=xpand_variables(x)
                                        
                    stud=regexp(x,['\<(',all_names,')\>'],'match');
                    
                    stud=[stud{:}];
                    
                    if ~isempty(stud)
                        
                        stud=stud{1};
                        
                        tmpx=x;
                        
                        x=cell(1,np);
                        
                        for ip=1:np
                            
                            x{ip}=regexprep(tmpx,stud,[stud,'_',panel_members{ip}]);
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
            
        function set_chain()
                        
            theMC=parser.initialize_markov_chain(mc_chains(ii).name,...
                mc_chains(ii).number_of_states,...
                'is_endogenous',~isempty(mc_chains(ii).endogenous_probabilities),...
                'param_list',cp(:).',...
                'is_switching',mc_chains(ii).number_of_states>1);
            
            if isempty(markov_chains)
                
                markov_chains=theMC;
                
            else
                
                markov_chains(1,ii)=theMC;
                
            end
            
        end
        
    end

end

function newc=expand_code(c)

newc={};

for ii=1:numel(c)
    
     newc=[newc,recursive_engine(c(ii))]; %#ok<AGROW>
    
    
end

    function ci=recursive_engine(ci)
    
    leftbrac=find(ci{1}=='[',1,'first');
    
    if isempty(leftbrac)
        
        return
        
    end
    
    rightbrac=find(ci{1}==']',1,'first');
    
    ll=ci{1}(1:leftbrac-1);
    
    rr=ci{1}(rightbrac+1:end);
    
    between=str2num(ci{1}(leftbrac:rightbrac)); %#ok<ST2NM>
    
    n=numel(between);
    
    cj=cell(1,n);
    
    ci=cell(1,0);
    
    for jj=1:n
        
        cj{jj}=sprintf('%s%0.0f%s',ll,between(jj),rr);
        
        ci=[ci,recursive_engine(cj(jj))]; %#ok<AGROW>
        
    end
    
        
    end

end

function x=decellize(x)

if iscell(x{1})
    
    x=[x{:}];
    
end

end

function x=recell2matize(x)

if isempty(x)
    
    x='';
    
else
    edn=parser.cell2matize(x);
    
    x=strrep(edn,'(','');
    
    x=strrep(x,')','');
    
end

end

function solution=do_solution(p,markov_chains,mc_detail,issvar,nvars,nx)

Qfunc=struct();

for ic=1:numel(markov_chains)
    
    Qfunc.(markov_chains(ic).name)=map_transition_matrix(markov_chains(ic));
    
end

solution=struct();

solution.var=@resolve_var;

solution.transition_matrix=@do_transition_matrix;

solution.full=@resolve;

    function [out,retcode]=resolve_var(M,~)
        
        retcode=0;
        
        nregs=size(M,2);
        
        out=struct();
        
        % useful for the reidentification of the SVAR
        out.nvars=nvars;
        
        out.nx=nx;
        
        if any(isnan(M(:)))
            
             retcode=2;
             
             return
             
        end
        
        for ireg=1:nregs
            
            varTerms=reshape(M(1:p.var_terms,ireg),nvars,[]);
            
            if issvar
                
                out.A(:,:,ireg)=varTerms;
                
                out.S0(:,:,ireg)=M(p.var_terms+(1:p.s_terms),ireg);
                
                [out.B(:,:,ireg),out.S(:,:,ireg)]=reduced_formize(...
                    out.A(:,:,ireg),out.S0(:,:,ireg),nx,nvars);
                
            else
                
                out.B(:,:,ireg)=varTerms;
                
                out.S(:,:,ireg)=ivech(M(p.var_terms+(1:p.s_terms),ireg));
                
            end
            
        end
        
    end
        
    function [out,retcode]=resolve(M,v)
        
        [out,retcode]=resolve_var(M,v);
        
        if ~retcode
            
            [out.Q,retcode]=do_transition_matrix(M,v);
            
        end
        
    end
        
    function [f,retcode]=do_transition_matrix(M,v)
        
        retcode=0;
        
        f=struct();
        
        Q=1;
        
        for ii=1:numel(markov_chains)
            
            tmp=Qfunc.(markov_chains(ii).name)(M(:,1),v);
            
            Q=kron(Q,tmp);
            
            if ~myisvalid(tmp)
                
                retcode=3;
                
                return
                
            end
            
            f.(markov_chains(ii).name)=tmp;
            
        end
        
        f.Q=Q;
        
        function ok=myisvalid(Q)
            
            ok = all(Q(:)>=0) && ...
                all(Q(:)<=1) && ...
                all(abs(sum(Q,2)-1)<1e-10);
            
        end
        
    end

    function out=map_transition_matrix(mc)
        
        cn=mc.name;
        
        tpv=@(varargin)[];
        
        if ~isempty(mc_detail)
            
            pos=find(strcmp(cn,{mc_detail.name}));
            
            if ~isempty(pos)
                
                tpn=mc_detail(pos).endogenous_probabilities{1};
                
                tpv=mc_detail(pos).endogenous_probabilities{2};
                
            end
            
        end
        
        nstates=mc.number_of_states;
        
        theMap=nan(nstates^2,3);
        
        iter=0;
        
        for istate=1:nstates
            
            for jstate=1:nstates
                
                if istate==jstate,continue,end
                
                iter=iter+1;
                
                tp=sprintf('%s_tp_%0.0f_%0.0f',cn,istate,jstate);
                
                tploc=find(strcmp(tp,tpn));
                
                theMap(iter,:)=[tploc,istate,jstate];
            
            end
            
        end
        
        theMap=abstvar.encode_map(theMap(1:iter,:),[nstates,nstates]);
        
        out=memoize_transition_matrix(theMap,nstates,tpv);
        
    end

end

function [B,S]=reduced_formize(A,S0,nx,nvars)

A0=A(:,nx+(1:nvars),:);

CA=A(:,[1:nx,nx+nvars+1:end],:);

iA0=A0\eye(nvars);

B=iA0*CA;

S=iA0*diag(S0);

% S is the covariance matrix
S=S*S.';

end

function out=memoize_transition_matrix(map,n,tpv)

out=@engine;

    function Q=engine(M,v)
        
        values=tpv(M(:,1),v);
                
        Q=zeros(n);
        
        if ~isempty(values)
            
            Q(map(:,2))=values(map(:,1));
            
        end
        
        diag_terms=(0:n-1)*n+(1:n);
        
        sq=sum(Q,2);
        
        Q(diag_terms)=1-sq;
        
    end

end

function vx=update_variable(v,v_names)

vx='\d+';

if ~isempty(v)
    
    if strcmp(v,':')
        
        % don't do anything
        
    elseif all(isstrprop(v,'digit'))
        
        vx=v;
        
    else
        
        vx=find(strcmp(v,v_names));
        
        if isempty(vx)
            
            error(['variable ',v,' not found among relevant variable names'])
            
        end
        
        vx=int2str(vx);
        
    end
    
end

end

function eqtnx=update_equation(eqtn)

eqtnx='\d+';

if ~isempty(eqtn)
    
    if strcmp(eqtn,':')
        
        % don't do anything
        
    else
        
        eqtnx=eqtn;
        
    end
    
end

end
