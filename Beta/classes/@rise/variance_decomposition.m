function [Vardec,obj,retcode]=variance_decomposition(obj,varargin)

% PURPOSE: Computes variance decompositions of a MSRE model
%-------------------------------------------------------------
% USAGE:
% where:
%
%
%
%
%-------------------------------------------------------------
% EXAMPLE:
%-------------------------------------------------------------
% LOG:
% 1.
%
%
%-------------------------------------------------------------
% SEE ALSO:
%-------------------------------------------------------------

% written [December 28, 2010] by:
% Junior Maih, Dept of Economics
% The Central Bank of Norway
% junior.maih@norges-bank.no
% this update [December 09, 2012]

Defaults=struct('vardec_ergodic',true,...
    'vardec_theoretical',true,...
    'vardec_periods',[1 4 8 16 40 100 200 400],...
    'vardec_shocks','');... % list of shocks to include: the others will have variance zero
    

if isempty(obj)
    Vardec=Defaults;
    return
end

nobj=numel(obj);
if nobj>1
    Vardec=cell(1,nobj);
    for iobj=1:nobj
        [Vardec{iobj},obj(iobj)]=variance_decomposition(obj(iobj),varargin{:});
    end
    return
end

obj=set_options(obj,varargin{:});

[obj,retcode]=solve(obj);	

% extract the state and transition matrices
if retcode
    error([mfilename,':: model could not be solved'])
end

oldoptions=obj.options;
thisDefault=Defaults;
Fields=fieldnames(thisDefault);
for ifield=1:numel(Fields)
    vi=Fields{ifield};
    if isfield(oldoptions,vi)
        thisDefault.(vi)=oldoptions.(vi);
    end
end
vardec_ergodic=thisDefault.vardec_ergodic;
vardec_theoretical  =thisDefault.vardec_theoretical  ;
vardec_periods  =thisDefault.vardec_periods  ;
k=max(100,max(vardec_periods));

% collect the state matrices
T=obj.T;
R=obj.R;
Q=obj.Q;
SS=obj.steady_state_and_balanced_growth_path(1:obj.NumberOfEndogenous(2),:);

[endo_nbr,exo_nbr,horizon,nstates]=size(R);
if vardec_ergodic
    [T,R,SS]=expected_state_matrices();
end

endo_names={obj.varendo.name};
exo_names={obj.varexo.name};

for ishock=1:numel(exo_names)
    shock=exo_names{ishock};
    if is_observed(shock)
        R(:,ishock,:,:)=0;
    end
end

Vardec=struct();
if vardec_theoretical
    if vardec_ergodic || nstates==1
        [Vinfi,Vi,retcode]=theoretical_vardec_engine(T,R);
        if ~retcode||nargout>2
        Vardec.conditional=vardec2rise_time_series(Vi);
        Vardec.infinity=vardec2rise_time_series(Vinfi);
        end
    else
        retcode0=0;
        for istate=1:nstates
            [Vinfi,Vi,retcode]=theoretical_vardec_engine(T(:,:,istate),R(:,:,:,istate));
            if ~retcode0
                retcode0=retcode;
            end
            if retcode && nargout<3
                break
            else
                Vardec(istate).conditional=vardec2rise_time_series(Vi);
                Vardec(istate).infinity=vardec2rise_time_series(Vinfi);
            end
        end
        retcode=retcode0;
    end
    if retcode && nargout<3
        error('variance could not be solved')
    end
else
    % compute from simulation
    error('Simulated variance decomposition not yet implemented')
end

    function db=vardec2rise_time_series(V)
        db=struct();
        this_k=size(V,3);
        for iendo=1:endo_nbr
            thisdec=transpose(reshape(V(iendo,:,:),exo_nbr,this_k));
            db.(endo_names{iendo})=rise_time_series(1,thisdec,exo_names);
        end
    end

    function flag=is_observed(shockname)
        flag=false;
        if ~isempty(obj.varobs_exo)
            obsNames={obj.varobs_exo.name};
            flag=ismember(shockname,obsNames);
        end
    end

    function [Vinfi,Vi,retcode0]=theoretical_vardec_engine(T,R)
        R=reshape(R,endo_nbr,exo_nbr*horizon);
        V=zeros(endo_nbr);
        Vi=zeros(endo_nbr,exo_nbr,k);
        Vii=zeros(endo_nbr,endo_nbr,exo_nbr);
        RR=R*R';
        for i=1:k
            V=T*V*T'+RR;
            [Vi(:,:,i),Vii]=decompose_variance(Vi(:,:,i),V,Vii);
        end
        Vinfi=zeros(endo_nbr,exo_nbr);
        [Vinf,retcode0]=lyapunov_equation(T,RR);
        if retcode0
            return
%             error('Variance could not be solved')
        end
        % deal with zero variances
        Vinfi=decompose_variance(Vinfi,Vinf);
        function [V,Vkk]=decompose_variance(V,total_variance,Vkk)
            total_variance=diag(total_variance);
            total_variance(total_variance<1e-12)=1;
            for iexo=1:exo_nbr
                Ri=zeros(size(R));
                locs=iexo:exo_nbr:exo_nbr*horizon;
                Ri(:,locs)=R(:,locs);
                RRi=Ri*Ri';
                if nargin<3
                    [V00,retcode1]=lyapunov_equation(T,RRi);
                    if retcode1
                        error('Variance could not be solved')
                    end
                    Vkk=V00;
                else
                    Vkk(:,:,iexo)=T*Vkk(:,:,iexo)*T'+RRi;
                    V00=Vkk(:,:,iexo);
                end
                V(:,iexo)=diag(V00)./total_variance;
            end
        end
    end

    function [A,B,SS_t]=expected_state_matrices()
        nreg=size(Q,1);
        probs_t=[eye(nreg)-Q';ones(1,nreg)]\[zeros(nreg,1);1];
        A=probs_t(1)*T(:,:,1);
        B=probs_t(1)*R(:,:,:,1);
        SS_t=probs_t(1)*SS(:,1);
        for st=2:nreg
            A=A+probs_t(st)*T(:,:,st);
            B=B+probs_t(st)*R(:,:,:,st);
            SS_t=SS_t+probs_t(st)*SS(:,st);
        end
    end

end
%%%%function Packeddec_db=SelectionAndRepackagingEngine(Vardec,Packages,var_list)
%%%%
%%%%endo_names=fieldnames(Vardec);
%%%%% indices of endogenous variables
%%%%i_var=locate_variables(var_list,char(endo_names));
%%%%out=setdiff(1:numel(endo_names),i_var);
%%%%Packeddec_db=rmfield(Vardec,endo_names(out));
%%%%
%%%%NumberOfPackages=size(Packages,2);
%%%%if NumberOfPackages==0
%%%%    return
%%%%end
%%%%%
%%%%ContributingNames=Packeddec_db.(deblank(var_list(1,:))).varnames;
%%%%startdate=Packeddec_db.(deblank(var_list(1,:))).start;
%%%%
%%%%% Now Repackage things up
%%%%DejaVu=[];
%%%%PackNames='';
%%%%Package_id=cell(1,NumberOfPackages);
%%%%for j=1:NumberOfPackages
%%%%    PackNames=strvcat(PackNames,Packages{1,j}); %#ok<VCAT>
%%%%    Package_id{j} = locate_variables(Packages{2,j},ContributingNames);
%%%%    tmp=intersect(DejaVu,Package_id{j});
%%%%    if ~isempty(tmp)
%%%%        tmp=mat2str(transpose(tmp(:)));
%%%%        error([mfilename,':: shocks ',tmp,' have been declared at least twice'])
%%%%    end
%%%%    DejaVu=union(DejaVu,Package_id{j});
%%%%end
%%%%Rest=setdiff(1:Packeddec_db.(deblank(var_list(1,:))).NumberOfVariables,DejaVu);
%%%%if ~isempty(Rest)
%%%%    PackNames=strvcat(PackNames,'Rest'); %#ok<VCAT>
%%%%    Package_id{NumberOfPackages+1}=Rest;
%%%%    NumberOfPackages=NumberOfPackages+1;
%%%%end
%%%%for j=1:numel(i_var)
%%%%    db=Packeddec_db.(deblank(var_list(j,:)));
%%%%    data=cell2mat(db.data(2:end,2:end));
%%%%    Newdata=nan(db.NumberOfObservations,NumberOfPackages);
%%%%    for p=1:NumberOfPackages
%%%%        Newdata(:,p)=sum(data(:,Package_id{p}),2);
%%%%    end
%%%%    Packeddec_db.(deblank(var_list(j,:)))=rise_time_series(startdate,Newdata,PackNames);
%%%%end
