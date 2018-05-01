function [ab,fval,retcode]=find_hyperparameters(space,cdfn,plb,pub,prob,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

number_of_starting_values=100;

debug=~true;

lb=space(:,1);

ub=space(:,2);

x0=initial_values();

happy_ending=@(exitflag,fval)ismember(exitflag,[1,3,5])&&fval<=1e-6;
% happy_endings=[1 % FMINCON: First order optimality conditions satisfied.|| LSQNONLIN: converged to a solution.
%     3]; % FMINCON: Change in objective function too small.|| LSQNONLIN : Change in RESNORM too small 
objective=@distributions.hyperparameter_residuals;

options=optimset('tolfun',1e-12,'tolx',1e-12);

if debug
    
    options.Display='iter';
    
else
    
    options.Display='off';%,'Algorithm','interior-point'
    
end

retcode=0;

done=false;

iter=0;

while ~done
    
    iter=iter+1;
    
    [ab,fval,~,exitflag]=lsqnonlin(@rich_man,x0(:,iter),lb,ub,options,'lsqnonlin');
    
    fval=norm(fval);
    
    he=happy_ending(exitflag,fval);
    
    if ~he
        
        [ab,fval,exitflag]=fsolve(@poor_man,x0(:,iter),options,'lsqnonlin');
        
        fval=norm(fval);
        
        he=happy_ending(exitflag,fval);
        
        if ~he
            
            [ab,fval,exitflag]=fmincon(@rich_man,x0(:,iter),[],[],[],[],lb,ub,[],options,'fmincon');
            
            he=happy_ending(exitflag,fval);
            
            if ~he
                
                [ab,fval,exitflag]=fminsearch(@poor_man,x0(:,iter),options,'fminsearch');
                
                he=happy_ending(exitflag,fval);
                
            end
            
        end
        
    end
    
    done=he || iter>=number_of_starting_values;
    
end

if ~he
    
    retcode=402;
    
end

if debug
    
    disp(ab)
    
    disp(fval)
    
    disp(retcode)
    
end

    function fval=rich_man(x,callingfunc)
        
        fval=objective(x,cdfn,plb,pub,prob,callingfunc,varargin{:});
        
    end

    function fval=poor_man(x,callingfunc)
        
        if any(x<lb)||any(x>ub)
            
            fval=1e+8;
            
        else
            
            fval=objective(x,cdfn,plb,pub,prob,callingfunc,varargin{:});
            
        end
        
    end

    function x0=initial_values()
        
        x0=[[plb,pub]',rand(size(lb,1),number_of_starting_values-1)];
        
        LB=lb(:,ones(1,number_of_starting_values));
        
        UB=ub(:,ones(1,number_of_starting_values));
        
        badneg=x0<LB;
        
        badpos=x0>UB;
        
        while any(any(badneg))||any(any(badpos))
            
            x0(badneg)=x0(badneg)+rand;
            
            x0(badpos)=x0(badpos)-rand;
            
            badneg=x0<LB;
            
            badpos=x0>UB;
            
        end
        
    end

end


