function [T,R,Z,H,Q,sstate,init,growth]=state_space_wrapper(syst)

h=numel(syst.T);

nshocks=size(syst.Te{1},2);

T=zeros(syst.m,syst.m,h);

R=zeros(syst.m,nshocks,h);

Z=syst.obs_id;

nobs=numel(Z);

H=zeros(nobs,nobs,h);

init=struct();
init.a=cell2mat(syst.a);
init.P=zeros(syst.m,syst.m,h);

Q=repmat(eye(nshocks),[1,1,h]);

for ireg=1:h
    
    T(:,syst.state_vars_location,ireg)=syst.Tx{ireg};
    
    R(:,:,ireg)=syst.Te{ireg};
    
    if ~isempty(syst.H{ireg})
        
        H(:,:,ireg)=syst.H{ireg};
        
    end
    
    init.P(:,:,ireg)=syst.P{ireg};
    
end

sstate=cell2mat(syst.steady_state);

growth=imag(cell2mat(syst.Tsig));

end