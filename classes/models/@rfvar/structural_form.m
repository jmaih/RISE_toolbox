function struct_forms=structural_form(obj,initcall)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

persistent oldobj endo_nbr nlags exo_nbr nx regimes_number B R Qzero ...
    zero_restr myperm_lags myperm_irfs sign_restr irf_sample_max Qsign %tags
if isempty(obj)
    struct_forms=struct('restrict_lags',{{}},...
		'restrict_irf_sign',{{}},'restrict_irf_zero',{{}},...
        'structural_shocks',{{}},'irf_sample_max',10000);
    return
end

if nargin<2
    initcall=true;
end

if isempty(obj.options.structural_shocks)
    struct_forms=obj;
    if obj.options.debug
        disp('structural shocks and restrictions not declared for identification')
    end
    return
end

if initcall
    myperm_lags=@(x)transpose(x);
    myperm_irfs=@(x)x;
    obj=check_identification(obj);
    if obj.identification.over
        error('error the model is over-identified and cannot be rotated')
    end
    
    endo_nbr=obj.endogenous.number(1);
    nlags=obj.nlags;
    nx=obj.nx;
    exo_nbr=obj.exogenous.number(1);
    regimes_number=obj.markov_chains.regimes_number;
    irf_sample_max=obj.options.irf_sample_max;
    
    % get solution
    %-------------
    [B,~,R]=vartools.resolve(obj.solution,nlags,regimes_number);
    
    % build the F matrix
    %-------------------
    rest_names={obj.nonlinear_restrictions.name};
    zero_restr=strcmp(rest_names,'lag_struct_Then_irf_zero_restr');
    sign_restr=obj.nonlinear_restrictions(~zero_restr);
    zero_restr=obj.nonlinear_restrictions(zero_restr);
    % Here we don't even need to sort this... usefull only for checking
    % identification
    Qzero=zero_restr.Q(1,:);% [Qzero,~,tags]=rfvar.sort_Q(zero_restr.Q);
    Qsign=sign_restr.Q(1,:);
    oldobj=obj;
end

% struct_forms=rfvar.empty(0,1); %struct('A0',{},'Aplus',{});
% One can either write a loop around this function or change the default
% number of rotations...
%-------------------------------------------------------------------------
success=false;
rounds=0;
while ~success && rounds<irf_sample_max
    % rotate R to get a new impact matrix C
    %--------------------------------------
    if oldobj.identification.exact
        C=R;
    else
        C=generate_new_impact_matrix();
    end
    
    Fzero=cell(1,regimes_number);
    is_lag_structure=true;
    for io=1:numel(zero_restr.orders)
        thisOrder=zero_restr.orders{io};
        for istate=1:regimes_number
            if is_lag_structure
                add_lag_structure_zero_restrictions()
            else
                add_irf_zero_restrictions()
            end
        end
        is_lag_structure=~is_lag_structure;
    end
    
    n=size(Fzero{1},2);
    Cnew=C;
    for istate=1:regimes_number
        P=var_rotation(Fzero{istate});
        % Compute corresponding short-run impact
        %---------------------------------------
        Cnew{istate} = C{istate}*P;
        Cnew{istate}(abs(Cnew{istate})<1e-12)=0;
    end
    
    success=isempty(sign_restr.f);
    if ~success
        % now check the sign restrictions
        %--------------------------------
        Fsign=cell(1,regimes_number);
        thisOrder=sign_restr.orders{1};
        inner_success=true;
        for istate=1:regimes_number
            if inner_success
                add_irf_sign_restrictions()
                % we cannot rotate the VAR, so we just check that the signs
                % check if not we go out
                %----------------------------------------------------------
                for iq=1:numel(Qsign)
                    inner_success=isempty(Qsign{iq})||...
                        all(Qsign{iq}*Fsign{istate}(:,iq)>0);
                    if ~inner_success
                        break
                    end
                end
            end
        end
        success=inner_success;
    end
    rounds=rounds+1;
end
if success
    fprintf(1,'successful rotation found after %0.0f rounds\n',rounds);
else
    error(['could not find suitable rotation after ',...
        sprintf('%0.0f',irf_sample_max),' iterations'])
end
% push the successful object
%---------------------------
for istate=1:regimes_number
    [oldobj.solution.omg{istate},oldobj.solution.sig{istate}]=vartools.decompose_impact(Cnew{istate});
end
struct_forms=oldobj;

    function C=generate_new_impact_matrix()
        C=cell(1,regimes_number);
        for st=1:regimes_number
            newmatrix = randn(exo_nbr);
            [QQ,RR] = qr(newmatrix);
            for ii = 1:exo_nbr
                if RR(ii,ii)<0
                    QQ(:,ii) = -QQ(:,ii);
                end
            end
            C{st} = transpose(R{st}'*QQ);
        end
    end
    function add_irf_sign_restrictions
        blocks=thisOrder;
        nrest=numel(blocks);
        if nrest
            L0=Cnew{istate};
            is_infinite=blocks(end)==inf;
            if is_infinite
                blocks(end)=[];
            end
            nrest=numel(blocks);
            if nrest
                max_order=thisOrder(end);
                if blocks(1)==0
                    Fsign{istate}=[Fsign{istate};myperm_irfs(L0)];
                end
                T=vartools.companion_form(B{istate},obj.nlags,endo_nbr,nx);
                y0=[L0;zeros(endo_nbr*(nlags-1),exo_nbr)];
                for iorder=1:max_order
                    y1=T*y0;
                    if ismember(iorder,blocks)
                        Li=y1(1:endo_nbr,:);
                        Fsign{istate}=[Fsign{istate};myperm_irfs(Li)];
                    end
                    y0=y1;
                end
            end
            if is_infinite
                Linf=vartools.find_long_run(B{istate},L0,obj.nlags,endo_nbr);
                Fsign{istate}=[Fsign{istate};myperm_irfs(Linf)];
            end
        end
    end
    function add_irf_zero_restrictions
        blocks=thisOrder;
        nrest=numel(blocks);
        if nrest
            L0=C{istate};
            is_infinite=blocks(end)==inf;
            if is_infinite
                blocks(end)=[];
            end
            nrest=numel(blocks);
            if nrest
                max_order=thisOrder(end);
                if blocks(1)==0
                    Fzero{istate}=[Fzero{istate};myperm_irfs(L0)];
                end
                T=vartools.companion_form(B{istate},obj.nlags,endo_nbr,nx);
                y0=[L0;zeros(endo_nbr*(nlags-1),exo_nbr)];
                for iorder=1:max_order
                    y1=T*y0;
                    if ismember(iorder,blocks)
                        Li=y1(1:endo_nbr,:);
                        Fzero{istate}=[Fzero{istate};myperm_irfs(Li)];
                    end
                    y0=y1;
                end
            end
            if is_infinite
                Linf=vartools.find_long_run(B{istate},L0,obj.nlags,endo_nbr);
                Fzero{istate}=[Fzero{istate};myperm_irfs(Linf)];
            end
        end
    end
    function add_lag_structure_zero_restrictions()
        nrest=numel(thisOrder);
        if nrest
            max_order=thisOrder(end);
            A0=C{istate}\eye(endo_nbr);
            if thisOrder(1)==0
                Fzero{istate}=[Fzero{istate};myperm_lags(A0)];
            end
            Aplus=A0*B{istate}(:,1:endo_nbr*max_order);
            for iorder=1:max_order
                if ismember(iorder,thisOrder)
                    Ai=Aplus(:,1:endo_nbr);
                    Fzero{istate}=[Fzero{istate};myperm_lags(Ai)];
                end
                Aplus=Aplus(:,endo_nbr+1:end);
            end
        end
    end
    function P=var_rotation(fa0ap)
        P=zeros(n);
        for jj=1:n
            Qj=Qzero{1,jj};% <-- Qj=Qzero{1,tags(jj)}; % we don't even need to sort anything here
            Qtilde=[Qj*fa0ap
                P(:,1:jj-1)'];% P'
            [q,~]=qr(Qtilde');
            P(:,jj)=q(:,end);
        end
    end
end


