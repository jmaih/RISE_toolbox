function [d,shortcut_d]=format_dirichlet(d,vals,estim_names,shortcut_d)

if nargin<4
    
    shortcut_d=[];
    
    if nargin<3
        
        estim_names=[];
        
        if nargin<2
            
            vals=[];
            
            if nargin<1
                
                d=[];
                
            end
            
        end
        
    end
    
end


if isempty(d)
    
    if ~isempty(shortcut_d)
        
        error('newd must be empty when d is empty')
        
    end
    
    [d,shortcut_d]=template();
    
    if isempty(vals)
        
        return
        
    end
    
end

nd=numel(d);

last_id=nd+1;

s02=vals{1}.^2;

vals=reshape(vals(2:end),2,[]);

pnames=vals(1,:);

m_main=cell2mat(vals(2,:));

m0=1-sum(m_main);

a_sum=m0*(1-m0)/s02-1;

if a_sum<=0
    
    error(['dirichlet # ',int2str(nd+1),...
        ' appears to have a too big standard deviation'])
    
end

m=[m_main(:);m0];

h=numel(m);

a=a_sum*m;

[d(last_id).a,d(last_id).b,d(last_id).moments,d(last_id).fval,...
    d(last_id).space]=distributions.dirichlet(a);

d(last_id).pointers=1:numel(pnames);

d(last_id).n_1=h-1;

d(last_id).names=pnames;

if ~isempty(estim_names)
    % where the estimates will be pushed
    loc=locate_variables(pnames,estim_names);
    
    d(last_id).location=loc;
    
    shortcut_d(last_id)=utils.distrib.dirichlet_shortcuts(d(last_id).a,...
        d(last_id).location,[],[]);
    
end

end

function [d,newd]=template()

d=struct('a',{},'b',{},'moments',{},'n_1',{},...
    'pointers',{},'location',{},'names',{});

newd=utils.distrib.dirichlet_shortcuts();

end