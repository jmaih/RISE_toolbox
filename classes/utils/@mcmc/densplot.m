function hdl=densplot(obj,pname,chain_id,N)

if nargin<4
    
    N=[];
    
    if nargin<3
        
        chain_id=[];
        
    end
    
end

if isempty(N)
    
    N=250;
    
end

x=load_draws(obj,pname,chain_id);

x=vec(x.');

[f_kdens,x_kdens]=distributions.kernel_density(x,[],[],'normal',N);

hdl0=plot(x_kdens,f_kdens);

title(pname)

axis tight

if nargout
    
    hdl=hdl0;
    
end

end