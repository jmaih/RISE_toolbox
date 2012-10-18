function eqn_obj=symbolic_derivator(endo_names,exo_names,par_names,lead_lag_incidence,equations,order)
% the dictionary contains all the names of the parameters, the endogenous,
% the exogenous

% 1- replace leads and lags by something else

if ischar(equations)
    equations=cellstr(equations);
end
number_of_equations=numel(equations);
% remove possible semi-colons in each equation
for eq=1:number_of_equations
    equations{eq}=strrep(equations{eq},';','');
end

if ischar(endo_names)
    endo_names=cellstr(endo_names);
end
if ischar(exo_names)
    exo_names=cellstr(exo_names);
end
if ischar(par_names)
    par_names=cellstr(par_names);
end
dictionary=union(union(endo_names,exo_names),par_names);
dictionary=dictionary(:);

% augment the dictionary with lags and leads using information in the
% lead_lag_incidence matrix. Replace (-1) with __L and (+1) with __F.
lagged_vars_id=find(lead_lag_incidence(:,1));
nlags=numel(lagged_vars_id);
auxiliary_lagged_vars=cell(nlags,1);
old_lagged_vars=cell(nlags,1);
for ii=1:nlags
    oldsubstr1=strcat(endo_names{lagged_vars_id(ii)},'(-1)');
    oldsubstr2=strcat(endo_names{lagged_vars_id(ii)},'{-1}');
    newsubstr=strcat(endo_names{lagged_vars_id(ii)},'__L');
    auxiliary_lagged_vars{ii}=newsubstr;
    old_lagged_vars{ii}=oldsubstr2;
    % replace the equations simultaneously
    for eq=1:number_of_equations
        equations{eq}=strrep(equations{eq},oldsubstr1,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr2,newsubstr);
    end
end
lead_vars_id=find(lead_lag_incidence(:,3));
nleads=numel(lead_vars_id);
auxiliary_lead_vars=cell(nleads,1);
old_lead_vars=cell(nleads,1);
for ii=1:nleads
    oldsubstr1=strcat(endo_names{lead_vars_id(ii)},'(+1)');
    oldsubstr2=strcat(endo_names{lead_vars_id(ii)},'{+1}');
    oldsubstr3=strcat(endo_names{lead_vars_id(ii)},'(1)');
    oldsubstr4=strcat(endo_names{lead_vars_id(ii)},'{1}');
    newsubstr=strcat(endo_names{lead_vars_id(ii)},'__F');
    auxiliary_lead_vars{ii}=newsubstr;
    old_lead_vars{ii}=oldsubstr4;
    % replace the equations simultaneously
    for eq=1:number_of_equations
        equations{eq}=strrep(equations{eq},oldsubstr1,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr2,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr3,newsubstr);
        equations{eq}=strrep(equations{eq},oldsubstr4,newsubstr);
    end
end
% add the auxiliary variables to the dictionary
dictionary=[auxiliary_lagged_vars;dictionary;auxiliary_lead_vars];
% list the variables to differentiate with respect to, ordering them
% according to lead_lag_incidence and adding exogenous at the end
with_respect_to=[auxiliary_lagged_vars
    endo_names(:)
    auxiliary_lead_vars
    exo_names(:)];

auxiliary_vars=[auxiliary_lagged_vars
    auxiliary_lead_vars];
old_vars=[old_lagged_vars
    old_lead_vars];

% create symbolic expressions for all the elements in the dictionary
for ii=1:numel(dictionary)
    eval([dictionary{ii},'=sym(dictionary{ii});'])
end

% symbolize every equation
for eq=1:number_of_equations
    equations{eq}=sym(equations{eq});
end

% create equation objects
number_of_variables=numel(with_respect_to);
the_zero=sym('0');

% % % % % % % % equations0=equations;
% % % % % % % % ii=0;
% % % % % % % % while ii<order
% % % % % % % %     ii=ii+1;
% % % % % % % %     eq_nbr=numel(equations0);
% % % % % % % %     str=['order_',int2str(ii)];
% % % % % % % %     Batch=cell(eq_nbr,number_of_variables);
% % % % % % % %     deriv.(str)=cell(eq_nbr,number_of_variables);
% % % % % % % %     eq=0;
% % % % % % % %     while eq<eq_nbr
% % % % % % % %         eq=eq+1;
% % % % % % % %         vv=0;
% % % % % % % %         while vv<number_of_variables
% % % % % % % %             vv=vv+1;
% % % % % % % %             if isequal(equations0{eq},the_zero)
% % % % % % % %                 tmp='0';
% % % % % % % %                 Batch{eq,vv}=the_zero;
% % % % % % % %             else
% % % % % % % %                 Batch{eq,vv}=diff(equations0{eq},with_respect_to{vv});
% % % % % % % %                 tmp=char(Batch{eq,vv});
% % % % % % % %                 for aa=1:nlags+nleads
% % % % % % % %                     tmp=strrep(tmp,auxiliary_vars{aa},old_vars{aa});
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %             deriv.(str){eq,vv}=tmp;
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % %     equations0=Batch(:);
% % % % % % % % end
%

% take one equation, find all the derivates and put them into an equation
% object
eq=0;
eq_nbr=numel(equations);
eqn_obj=rise_equation.empty(0);
 x={'first','second','third'};
while eq<eq_nbr
    eq=eq+1;
    ii=0;
    equations0=equations{eq};
    eqn_obj(eq,1)=rise_equation(equations0,with_respect_to,par_names,eq);
    while ii<order
        ii=ii+1;
        inner_eq_nbr=numel(equations0);
        i_eq=0;
        Batch=cell(inner_eq_nbr,number_of_variables);
        EQ_eq_ORDER_ii=cell(inner_eq_nbr,number_of_variables);
        while i_eq<inner_eq_nbr
            i_eq=i_eq+1;
            is_zero=isequal(equations0{i_eq},the_zero);
            vv=0;
            tmp='0';
            while vv<number_of_variables
                vv=vv+1;
                if is_zero
                    Batch{i_eq,vv}=the_zero;
                else
                    Batch{i_eq,vv}=diff(equations0{i_eq},with_respect_to{vv});
                    tmp=char(Batch{i_eq,vv});
                    for aa=1:nlags+nleads
                        tmp=strrep(tmp,auxiliary_vars{aa},old_vars{aa});
                    end
                end
                EQ_eq_ORDER_ii{i_eq,vv}=tmp;
            end
        end
        prop=[x{ii},'_order_derivatives'];
        eqn_obj(eq,1)=eqn_obj(eq,1).set_properties(prop,EQ_eq_ORDER_ii);
        equations0=Batch(:);
    end
end
% equation
% expression
% first_order
% second_order

