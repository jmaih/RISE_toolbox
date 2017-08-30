function [fc,viol,msg]=eval_objective_and_restrictions(c,objective,...
    restrictions,opt,varargin)

violObjective=[];

rcode=0;

if opt.restrictions_in_objective
    
    [fc,rcode,violObjective]=objective(c,varargin{:});
    
elseif opt.returns_retcode
    
    [fc,rcode]=objective(c,varargin{:});
    
else
    
    fc=objective(c,varargin{:});
    
end

msg=utils.error.decipher(rcode);

viol1=restrictions(c);

viol=[violObjective;viol1];

viol=utils.estim.penalize_violations(viol,[],opt.restrictions_same_weights);

end