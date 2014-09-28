function [x,f,viol,funevals]=generate_candidates(objective,lb,ub,n,restrictions,penalty,varargin)
if isempty(restrictions)
    restrictions=@(z)[];
end
npar=size(lb,1);
x=nan(npar,n);
f=nan(1,n);
fcount=nan(1,n);
viol=cell(1,n);
MaxIter=50;
max_trials=50;
funevals=0;
% can we get 4 outputs?
success=nargout(objective)>=2;
msg='';
the_loop=@loop_body;
for ii=1:n
    [x(:,ii),f(ii),viol{ii},fcount(ii)]=the_loop(ii);
end
funevals=funevals+sum(fcount);

    function [xi,fi,violi,fcount]=loop_body(ii)
        invalid=true;
        iter=0;
        fcount=0;
        while invalid
            if iter>=MaxIter
                error([mfilename,':: could not generate a valid candidate after ',...
                    int2str(MaxIter*max_trials),' attempts'])
            end
            [xi,fi,violi]=draw_and_evaluate_vector();
            iter2=0;
            while iter2<max_trials && sum(violi)>0
                iter2=iter2+1;
                [c,fc,violc]=draw_and_evaluate_vector();
                if sum(violc)<sum(violi)
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
        viol=restrictions(c);
        viol=viol(viol>0);
        if isempty(viol)
            viol=0;
        end
    end

end