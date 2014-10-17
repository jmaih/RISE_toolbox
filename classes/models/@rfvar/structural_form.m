function newobj=structural_form(obj,varargin)
% structural_form finds A structural form given the imposed restrictions
%
% Syntax
% -------
% ::
%
%   newobj=structural_form(obj)
%   newobj=structural_form(obj,varargin)
%
% Inputs
% -------
%
% - **obj** : [rfvar] : reduced form VAR object
%
% - varargin : standard optional inputs **coming in pairs**. Among which:
%
%   - **restrict_lags** \: [cell array|{''}] \: restrictions on the lag
%       structure. There are two equivalent syntaxes for this:
%
%       - {'var_name1@var_name2{lag}'}
%       - {'alag(var_name1,var_name2)'} \: here alag should be understood as
%         a-lag, where lag is the "lag" e.g. a1(infl,unemp) means unemp
%         does not enter the infl equation at lag 1.
%   - **restrict_irf_sign** \: [cell array|{''}] \: sign restrictions on the
%       impulse responses. The general syntax is
%       {'var_name{period}@shock_name','sign'} and the default period is
%       "0" (for contemporaneous). That means
%       {'var_name{0}@shock_name','+'} and {'var_name@shock_name','+'}
%       are equivalent
%   - **restrict_irf_zero** \: [cell array|{''}] \: zero restrictions on the
%       impulse responses. The general syntax is
%       {'var_name{period}@shock_name'} and the default period is
%       "0" (for contemporaneous). That means
%       {'var_name{0}@shock_name'} and {'var_name@shock_name'}
%       are equivalent
%   - **structural_shocks** \: [cell array|{''}] \: List of structural
%       shocks. The shock names can be entered with or without their
%       description. For instance : 
%       - {'E_PAI','E_U','E_MP'}
%       - {'E_PAI','"inflation shock"','E_U','"unempl shock"','E_MP'}
%   - **irf_sample_max** : [numeric|{10000}] : maximum number of trials in
%       the drawing of rotation matrices
%
% Outputs
% --------
%
% - newobj : [rfvar]: new rfvar object with the drawn structural form
%
% More About
% ------------
%
% - RISE automatically orders the endogenous variables alphabetically and
%   tags each equation with one of the endogenous variables. This may be
%   useful for understanding the behavior of **restrict_lags** above.
%
% - The Choleski identification scheme is not implemented per se. The user
%   has to explicitly enter the zeros in the right places. This gives the
%   flexibility in implementing the restrictions. For instance, one could
%   imagine a scheme in which choleski restrictions hold only in the long
%   run.
%
% - With only zero restrictions, one cannot expect the impulse responses to
%   automatically have the correct sign. The rotation imposes zero
%   restrictions but not the sign. If you would like to have
%   correctly-signed impulse responses there are two choices:
%   - explicitly add sign restrictions
%   - multiply the impulse responses for the wrongly-signed shock with
%   minus.
%
% - Many periods can be entered simultaneously. For instance 
%   'var_name{0,3,5,10:20,inf}@shock_name' 
%
% - long-run restrictions are denoted by "inf". For instance 
%   'var_name{inf}@shock_name' 
%
% - Identification for Markov switching VARs is not implemented/supported. 
%
% Examples
% ---------
%
% See also:

if isempty(obj)
    newobj=struct('restrict_lags',{{}},...
        'restrict_irf_sign',{{}},'restrict_irf_zero',{{}},...
        'structural_shocks',{{}},'irf_sample_max',10000);
    return
end

if isempty(obj.options.structural_shocks)
    newobj=obj;
    if obj.options.debug
        disp('structural shocks and restrictions not declared for identification')
    end
    return
end

obj=set(obj,varargin{:});

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

if regimes_number>1
    error('identification for Markov-switching VARs is not implemented')
end

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
Qsign=sign_restr.Q(1,:);
[Qzero,~,tags_zero]=rfvar.sort_Q(zero_restr.Q);
itags_zero(tags_zero)=1:numel(tags_zero);

% newobj=rfvar.empty(0,1); %struct('A0',{},'Aplus',{});
% One can either write a loop around this function or change the default
% number of rotations...
%-------------------------------------------------------------------------
success=false;
rounds=0;
while ~success && rounds<irf_sample_max
    % rotate R to get a new impact matrix C
    %--------------------------------------
    if obj.identification.exact
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
    
    number_of_zero_restrictions=size(Fzero{1},2);
    Cnew=C;
    for istate=1:regimes_number
        P=var_rotation(Fzero{istate});
        % Compute corresponding short-run impact
        %---------------------------------------
        Cnew{istate} = C{istate}*P(:,itags_zero);
        Cnew{istate}(abs(Cnew{istate})<1e-12)=0;
        if obj.options.debug
            check_rotation_adequacy(C{istate},Cnew{istate});
            bigtable=cell(obj.endogenous.number(end)+1,sum(obj.exogenous.number)+1);
            bigtable(2:end,1)=obj.endogenous.name;
            bigtable(1,2:end)=obj.exogenous.name;
            bigtable(2:end,2:end)=num2cell(Cnew{istate});
            disp(bigtable)
            keyboard
        end
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
    [obj.solution.omg{istate},obj.solution.sig{istate}]=vartools.decompose_impact(Cnew{istate});
end
newobj=obj;

    function check_rotation_adequacy(C,Cnew)
        SIG=C*C';
        SIG_new=Cnew*Cnew';
        discrepancy=max(abs(SIG(:)-SIG_new(:)));
        if discrepancy>1e-9
            error('rotation code does not work')
        end
    end

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
        if number_of_zero_restrictions
            P=zeros(number_of_zero_restrictions);
            for jj=1:number_of_zero_restrictions
                Qj=Qzero{1,jj};
                Qtilde=[Qj*fa0ap
                    P(:,1:jj-1)'];
                [q,~]=qr(Qtilde');
                P(:,jj)=q(:,end);
            end
        else
            P=eye(endo_nbr);
        end
    end
end


