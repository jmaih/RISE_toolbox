function db=growth_database(obj,endhist_date,end_sample,growth_type)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)

    db=struct();

    return

end

if ~ismember(growth_type,{'zero','steady'})
    
    error('growth_type should be "zero" or "steady"')
    
end

zero_growth_=strcmp(growth_type,'zero');

% needs to be adapted for VARs
sv=get(obj,'state_vars');

maxlag=max(cell2mat(sv(2,:)))+1;

is_log_var=obj.endogenous.is_log_var;

ss=get(obj,'sstate');

g=get(obj,'growth');

endhist_date=date2serial(endhist_date);

starthist_date=endhist_date-maxlag+1;

last_date=date2serial(end_sample);

nperiods=numel(starthist_date:last_date);

db=struct();

nv=obj.endogenous.number;

span=(1-maxlag:nperiods-maxlag)';

h=obj.markov_chains.regimes_number;

d=nan(numel(span),h);

if zero_growth_
    span(:)=0;
end

for iv=1:nv
    
    name=obj.endogenous.name{iv};
    
    do_growth();
            
    db.(name)=ts(starthist_date,d);
    
end

d=zeros(nperiods,h);

for ix=1:sum(obj.exogenous.number)
    
    name=obj.exogenous.name{ix};
    
    db.(name)=ts(starthist_date,d);
    
end

    function do_growth()
        
        ssv=ss.(name);
        
        gv=g.(name);
        
        for ireg=1:h
            
            if is_log_var(iv)
                
                d(:,ireg)=ssv(ireg)*gv(ireg).^span;
                
            else
                
                d(:,ireg)=ssv(ireg)+gv(ireg)*span;
                
            end
            
        end 
        
    end

end