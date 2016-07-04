function hdl=scatterplot(obj,pname1,pname2,chain_id,varargin)

if nargin<4
    
    chain_id=[];
    
end

x1=load_draws(obj,pname1,chain_id);

x1=vec(x1.');

x2=load_draws(obj,pname2,chain_id);

x2=vec(x2.');

hdl0=scatter(x1,x2,varargin{:});

xlabel(pname1)

ylabel(pname2)

axis tight

if nargout
    
    hdl=hdl0;
    
end

end