function obj=translate_restrictions(obj)
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

if isempty(obj)
    obj=struct();
    return
end
restr=struct('name',{},'orders',{},'array',{});
types=regexp(fieldnames(obj.options),'restrict_\w+','match');
types=[types{:}];
endo_names=obj.endogenous.name;
endo_nbr=numel(endo_names);
struct_shocks=get(obj,'exo_list(~observed)');
exo_nbr=obj.exogenous.number(1); % not observed!!!

for itype=1:numel(types)
    restr_name=types{itype};
    is_lag_structure=strcmp(restr_name,'restrict_lags');
    is_sign_restriction=strcmp(restr_name,'restrict_irf_sign');
    thisRestriction=obj.options.(restr_name);
    orders=[];
    array=nan(endo_nbr,exo_nbr,0);
    % if there are complicated lag structure restrictions they will
    % crash as there is no theory for checking the identification
    % implied by such restrictions. Hence if the model is a
    % reduced-form VAR I should check that complicated restrictions
    % involving operations are not provided.
    nrows=size(thisRestriction,1);
    if ~isempty(thisRestriction)
        % standardize all restrictions to Y@R{}
        %--------------------------------------
        % turn a3(infl,ouput) into infl@ouput{3}
        thisRestriction(:,1)=regexprep(thisRestriction(:,1),'a(\d+)\((\w+),(\w+)\)','$2@$3{$1}');
        
        left=thisRestriction(:,1);
        get_orders();
        no=numel(orders);
        array=nan(endo_nbr,exo_nbr,no);
        if is_lag_structure
            % get the template for the SVAR and push the block
            % exogenous restrictions
            %-------------------------------------------------
            svar_param_template=vartools.set_structural_matrices_structure('svar',...
                obj.endogenous.number,obj.nlags,obj.nx,obj.endogenous.is_block_exogenous);
            lag_names=regexp(svar_param_template(1,:),'(?<!\w+)a\d+(?!\w+)','match');
            lag_names=[lag_names{:}];
            lag_pos=locate_variables(lag_names,svar_param_template(1,:));
            lags=strrep(lag_names,'a','');
            lags=cell2mat(strcat(lags,','));
            lags=eval(['[',lags(1:end-1),']']);
            for io=1:no
                array(:,:,io)=svar_param_template{2,lag_pos(orders(io)==lags)};
            end
        end
        if is_sign_restriction
            right=thisRestriction(:,2);
        end
        fill_array();
    end
    restr(itype)=struct('name',restr_name,'orders',orders,'array',array);
end

% zero restrictions: mix lag structure restrictions with irf restrictions
%------------------------------------------------------------------------
names={restr.name};
fields={'restrict_lags','restrict_irf_zero'};
nf=numel(fields);
orders_=cell(1,nf);
f0=zeros(0,exo_nbr);
for iname=1:nf
    loc= strcmp(fields{iname},names);
    % permute columns and rows to be consistent with RWZ's notation
    %--------------------------------------------------------------
    is_irf=strcmp(fields{iname},'restrict_irf_zero');
    for ir=1:size(restr(loc).array,3)
        item=restr(loc).array(:,:,ir);
        if ~is_irf
            % in RWZ the irfs are transposed in my notation, it should
            % be the opposite
            item=transpose(item);
        end
        f0=[f0;item]; %#ok<*AGROW>
    end
    orders_{iname}=restr(loc).orders;
end
[Q0,f0]=find_Q(f0,0);
nonlinear_restrictions=struct('name','lag_struct_Then_irf_zero_restr',...
    'Q',{Q0},'orders',{orders_},'f',f0);

% sign restrictions
%------------------
loc=strcmp('restrict_irf_sign',names);
fs=zeros(0,exo_nbr);
for ir=1:size(restr(loc).array,3)
    % in RWZ the irfs are transposed in my notation, it should
    % be the opposite
    item=restr(loc).array(:,:,ir);
    fs=[fs;item];
end
orders_={restr(loc).orders};
[Qs,fs]=find_Q(fs,[1,-1]);
nonlinear_restrictions(2)=struct('name','sign_restr','Q',{Qs},...
    'orders',{orders_},'f',fs);

obj.nonlinear_restrictions=nonlinear_restrictions;

    function [Q,f]=find_Q(f,arg)
        nrf=size(f,1);
        zero_restr=numel(arg)==1;
        % for each shock build a Q matrix
        %--------------------------------
        Q=cell(2,exo_nbr);
        for ishock=1:exo_nbr
            if zero_restr
                loc=find(f(:,ishock)==arg);
            else
                loc=find(f(:,ishock)==arg(1)|f(:,ishock)==arg(2));
            end
            nloc=numel(loc);
            tmp=zeros(nrf,nloc);
            if zero_restr
                tmp(loc,:)=eye(nloc);
            else
                tmp(loc,:)=diag(f(loc,ishock));
            end
            Q{1,ishock}=transpose(tmp);
            if isempty(tmp)
                Q{2,ishock}=0;
            else
                Q{2,ishock}=rank(tmp);
            end
        end
        % final transformation on f:
        %---------------------------
        if zero_restr
            % for zero restrictions, set all nans to 1, so that the
            % conditions can be checked
            %--------------------------------------------------------------
            f(isnan(f))=1;
        else
            % for sign restrictions, set all nans to 0.5, so that the
            % conditions can be checked
            %--------------------------------------------------------------
            f(isnan(f))=0.5;
        end
    end
    function fill_array()
        for irow=1:nrows
            [oo,v,x]=do_one(left{irow});
            vloc=find(strcmp(v,endo_names));
            if isempty(vloc)
                error(['left-hand side variable "',v,'" could not be located'])
            end
            if is_lag_structure
                % I still don't get how the shocks will be mixed with the
                % lags...
                xloc=find(strcmp(x,endo_names));
            else
                xloc=find(strcmp(x,struct_shocks));
            end
            if isempty(xloc)
                error(['right-hand side variable "',x,'" could not be located'])
            end
            for io_=1:numel(oo)
                pos=find(oo(io_)==orders);
                if is_sign_restriction
                    splus=strcmp(right{irow},'+');
                    array(vloc,xloc,pos)=splus-~splus;
                else
                    array(vloc,xloc,pos)=0;
                end
            end
        end
    end
    function get_orders()
        for irow=1:nrows
            oo=do_one(left{irow});
            orders=[orders,oo(:)'];
        end
        orders=unique(orders);
    end
    function [oo,v,x]=do_one(item)
        oo=0;
        lcb=find(item=='{');
        rtype='}';
        if isempty(lcb)
            lcb=find(item=='(');
            rtype=')';
        end
        if ~isempty(lcb)
            rcb=find(item==rtype);
            if isempty(rcb)
                error('closing curly brace of parenthesis missing')
            end
            oo=item(lcb+1:rcb-1);
            if isempty(oo)
                error('order missing')
            end
            oo=sort(abs(eval(['[',oo,']'])));
        end
        if ~isequal(oo,unique(oo))
            error('redundant/repeated orders')
        end
        if nargout>1
            a=find(item=='@');
            if isempty(a)||numel(a)>1
                error(['unable to parse restriction "',item,'"'])
            end
            if isempty(lcb)
                v=item(1:a-1);
            else
                % remove the time as it can appear either before or after
                % the @ sign
                %---------------------------------------------------------
                item=[item(1:lcb-1),item(rcb+1:end)];
                a=find(item=='@');
                v=item(1:a-1);
            end
            x=item(a+1:end);
        end
    end
end