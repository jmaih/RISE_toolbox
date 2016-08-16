function x=one_step_simulator(sol,x0,e)

nstates=size(sol.Tx,3);

x=cell(1,nstates);

nrounds=size(e,3);

for ii=1:nstates
    
    x{ii}=one_simulation(sol.Tx(:,:,ii),sol.ss{ii},sol.Tsig(:,1,ii),sol.Te(:,:,ii));
    
end

    function xx=one_simulation(Tx,xbar,Tsig,Te)
        
        for iround=1:nrounds
            
            x1=xbar+Tx*(x0-xbar)+Tsig+Te*e(:,:,iround);
            
            if iround==1
                
                xx=x1(:,ones(1,nrounds));
                
            else
                
                xx(:,iround)=x1;
                
            end
            
        end
        
    end

end