function [x,f,retcode]=generate_starting_point(objective,lb,ub,restrictions)
if nargin<4
    restrictions=[];
end
npar=size(lb,1);
c=lb+(ub-lb).*rand(npar,1);
if ~isempty(restrictions)
    c=reprocess_parameter_vector(c);
end
[fc,retcode]=objective(c);
x=c;
f=fc;

    function [c,modified]=reprocess_parameter_vector(c)
        % this function attempts to reduce the number of rejection of randomly
        % drawn vectors in the presence of linear or nonlinear restrictions. When
        % there are many restrictions, randomly drawing a vector of parameters that
        % satisfies them all can be challenging.
        % INPUTS:
        %       restrictions is a n x 2 cell array. For each row, the first column is the
        % list of the parameters involved in the restrictions and the second column
        % is a function that checks that the restriction holds.
        %       lb: lower bound
        %       ub: upper bound
        %       c: the parameter vector to modify
        % If a parameter is involved in several restrictions, the reprocessing may
        % fail. If a restriction function admits more than just the vector of
        % parameters to check, one can construct a function handle such as
        % @(x)f(x,a), where a has been previously defined.
        % this function is called by generate_candidates
        modifs=0;
        for irest=1:size(restrictions,1)
            involved=restrictions{irest,1};
            involved=involved(:);
            the_func=restrictions{irest,2};
            done=the_func(c);
            n=numel(involved);
            iter=0;
            while ~done
                iter=iter+1;
                c(involved)=lb(involved)+rand(n,1).*(ub(involved)-lb(involved));
                done=the_func(c)||iter>500;
            end
            modifs=modifs+iter;
        end
        modified=logical(modifs);
    end
end