function x=one_step_simulator(sol,x0,e,do_parallel)

if nargin < 4
    
    do_parallel = false;
    
end

nstates=size(sol.Tx,3);

x=cell(1,nstates);

nrounds=size(e,3);

if do_parallel
    
    Tx=sol.Tx;
    
    ss=sol.ss;
    
    Tsig=sol.Tsig;
    
    Te=sol.Te;
    
    theSimulator=@one_simulation;
    
    parfor (ii=1:nstates,utils.parallel.get_number_of_workers)
        
        x{ii}=theSimulator(Tx(:,:,ii),ss{ii},Tsig(:,1,ii),Te(:,:,ii));
        
    end
    
else
    
    for ii=1:nstates
        
        x{ii}=one_simulation(sol.Tx(:,:,ii),sol.ss{ii},sol.Tsig(:,1,ii),sol.Te(:,:,ii));
        
    end
    
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