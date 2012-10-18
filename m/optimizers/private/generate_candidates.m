function [x,f,funevals]=generate_candidates(objective,lb,ub,n,restrictions,penalty,varargin)
npar=size(lb,1);
x=nan(npar,n);
f=nan(1,n);
MaxIter=5000;
funevals=0;
% can we get 4 outputs?
success=false;
try %#ok<TRYNC>
    [junk,junk,junk,o4]=objective(lb+(ub-lb).*rand(npar,1),varargin{:}); %#ok<NASGU>
    success=true;
end
msg='';
for ii=1:n
    invalid=true;
    iter=0;
    while invalid
        if iter>=MaxIter
            error([mfilename,':: could not generate a valid candidate after ',...
                int2str(MaxIter),' attempts'])
        end
        c=lb+(ub-lb).*rand(npar,1);
        if ~isempty(restrictions)
            c=reprocess_parameter_vector(c);
        end
        if success
            [fc,junk,junk,o4]=objective(c,varargin{:});
            msg=decipher_error(o4);
        else
            fc=objective(c,varargin{:});
        end
        funevals=funevals+1;
        if fc<penalty
            invalid=false;
        else
            fprintf(1,'%5s %3.0d/%3.0d %8s %5.0d %8s %s \n',...
                'warrior #',ii,n,'iter',iter,'pb',msg);
        end
        iter=iter+1;
    end
    x(:,ii)=c;
    f(ii)=fc;
end

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
            n_i=numel(involved);
            iter_i=0;
            while ~done
                iter_i=iter_i+1;
                c(involved)=lb(involved)+rand(n_i,1).*(ub(involved)-lb(involved));
                done=the_func(c)||iter_i>500;
            end
            modifs=modifs+iter_i;
        end
        modified=logical(modifs);
    end
end