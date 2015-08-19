function [x,f,viol,funevals]=generate_candidates(objective,lb,ub,n,...
    restrictions,penalty,varargin)
% generate_candidates -- generate candidates for optimization
%
% Syntax
% -------
% ::
%
%   [x,f,viol,funevals]=generate_candidates(objective,lb,ub,n,...
%       restrictions,penalty,varargin)
%   
% Inputs
% -------
%
% - **objective** [function_handle]: objective to minimize
%
% - **lb** [vector]: lower bound of the search space
%
% - **ub** [vector]: upper bound of the search space
%
% - **n** [integer]: number of candidates to generate
%
% - **restrictions** [empty|function_handle]: function evaluating the
% violations
%
% - **penalty** [numeric]: value functions exceeding this value in
% absolute value are assigned this value 
%
% - **varargin** []: additional input arguments for the objective function
%
% Outputs
% --------
%
% - **x** [d x n matrix]: parameter vectors
%
% - **f** [row vector]: value function for each parameter vector
%
% - **viol** [row vector]: violation penalty for each parameter vector
%
% - **funevals** [integer]: number of function evaluations
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(restrictions)
    restrictions=@(z)[];
end
npar=size(lb,1);
x=nan(npar,n);
f=nan(1,n);
fcount=nan(1,n);
viol=zeros(1,n);
MaxIter=50;
max_trials=50;
funevals=0;
% can we get 2 outputs?
success=false;
try %#ok<TRYNC>
    c=lb+(ub-lb).*rand(npar,1);
    [~,junk]=objective(c,varargin{:}); %#ok<NASGU>
    success=true;
end
% success=nargout(objective)>=2;
msg='';
the_loop=@loop_body;
for ii=1:n
    [x(:,ii),f(ii),viol(ii),fcount(ii)]=the_loop(ii);
end
funevals=funevals+sum(fcount);

    function [xi,fi,violi,fcount]=loop_body(ii)
        invalid=true;
        iter=0;
        fcount=0;
        while invalid
            if iter>=MaxIter
                error([mfilename,':: could not generate a valid candidate after ',...
                    int2str(MaxIter*max_trials),' attempts. Exiting with the best vector'])
            end
            [xi,fi,violi]=draw_and_evaluate_vector();
            iter2=0;
            while iter2<max_trials && violi>0
                iter2=iter2+1;
                [c,fc,violc]=draw_and_evaluate_vector();
                if violc<violi
                    xi=c;
                    fi=fc;
                    violi=violc;
                end
            end
            fcount=fcount+1;
            if fi<penalty
                invalid=false;
            else
                fprintf(1,'%5s %3.0d/%3.0d %8s %5.0d %8s %s \n',...
                    'warrior #',ii,n,'iter',iter,'pb',msg);
            end
            iter=iter+1;
        end
    end

    function [c,fc,viol]=draw_and_evaluate_vector()
        c=lb+(ub-lb).*rand(npar,1);
        if success
            [fc,o4]=objective(c,varargin{:});
            msg=utils.error.decipher(o4);
        else
            fc=objective(c,varargin{:});
        end
        viol=utils.estim.penalize_violations(restrictions(c));
    end

end