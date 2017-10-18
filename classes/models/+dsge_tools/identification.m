function J=identification(m,p,obslist,ar)

if nargin<4
    
    ar=[];
    
    if nargin<3
        
        obslist=[];
        
        if nargin<2
            
            p=[];
            
        end
        
    end
    
end

if isempty(ar),ar=0;end

if isempty(p),p=get(m,'mode');end

if isempty(obslist)
    
    obslist=m.observables.name;
    
    obspos=locate_variables(obslist,m.endogenous.name);
    
else
    
    obspos=m.observables.state_id;
    
end

plist=fieldnames(p);

n=numel(plist);

p0=load_start();

d=ar+1;

m=set(m,'autocov_ar',ar);

J=utils.numdiff.jacobian(@obsjac,p0);

    function [my,retcode]=obsjac(params)
        
        if size(params,2)>1
            
            error('parameters should be a colum vector')
            
        end
        
        p1=set_params(params);
        
        [obj,retcode]=solve(m,'parameters',p1);
        
        ss=obj.solution.ss{1};
        
        AC=theoretical_autocovariances(obj);
        
        mu=ss(obspos);
        
        ACy=vech_ize(AC);
        
        my=[mu.',ACy.'].';
        
    end

    function v=vech_ize(C)
        
        v=cell(1,d);
        
        for ii=1:d
            
            v{ii}=vech(C(obspos,obspos,ii));
            
        end
        
        v=cell2mat(v);
        
        v=v(:);
        
    end

    function p1=set_params(x)
        
        p1=p;
        
        for ii=1:n
            
            p1.(plist{ii})=x(ii);
            
        end
        
    end

    function p0=load_start()
        
        p0=nan(n,1);
        
        for ii=1:n
            
            p0(ii)=p.(plist{ii});
            
        end
        
    end

end