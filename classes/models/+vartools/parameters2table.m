function tab=parameters2table(A)
% INTERNAL FUNCTION
%

[nvars,acols]=size(A);

nlags=floor(acols/nvars);

nx=acols-nlags*nvars;

colseps=nvars*ones(1,nlags);

if nx

    colseps=[nx,colseps];

end

t=mat2cell(A,nvars,colseps);

create_names=@(a,b)utils.char.create_names([],a,b);

rownames=create_names('eq',nvars);

colnames=create_names('lag',nlags);

if nx

    colnames=['deterministic',colnames];

end

try

    tab=table(t{:},'VariableNames',colnames,'RowNames',rownames);

catch me

    disp(me)

    tab=[];

    for ii=1:numel(colnames)

        disp(['###########',colnames{ii},'###########'])

        disp([rownames(:),num2cell(t{ii})])

    end

end


end
