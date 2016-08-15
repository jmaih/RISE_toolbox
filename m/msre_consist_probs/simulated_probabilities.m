function c=simulated_probabilities(solution,x0,inner_probabilities,control_shocks)

nchk=size(inner_probabilities,1);

c=zeros(1,nchk);

if nchk == 0
    
    return
    
end

x=one_step_simulator(solution,x0,control_shocks);

for ichk = 1:nchk
    
    locs = inner_probabilities{ichk,2}.chain_loc;
    
    xx = cell2mat(x(locs));
    
    c(ichk)=process(inner_probabilities{ichk,2}.func);
    
end

    function prob=process(func)
        
        good=func(xx);
        
        prob=sum(good)/numel(good);
        
    end


end