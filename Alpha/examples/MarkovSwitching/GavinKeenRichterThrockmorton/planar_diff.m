clc
nwrt=5;
order=10;
Prev={};
incmnt=500;
for io=1:order
    nprev=max(1,size(Prev,2));
    Current=cell(2,incmnt); % 1- indexes, 2-pointers
    cur_size=incmnt;
    iter=0;
    for iprev=1:nprev
        if isempty(Prev)
            last_ind=[];
            wrt=1:nwrt;
            Pointer=1; % column to differentiate
        else
            last_ind=Prev{1,iprev};
            wrt=last_ind(end):nwrt;
            Pointer=Prev{2,iprev};
        end
        theLegend=false(1,nwrt);
        for ieq=1:neqtns
            thisLegend=theLegend;
            starters=ispresent(wrt);
            thisLegend(starters)=true;
            if ispresent(Pointer)
                d=diff(eqtn(ieq,ibatch),starters,Pointer);
            end
        end
        Next_Pointers=1:numel(wrt);
        for irt=Next_Pointers
            iter=iter+1;
            if iter==cur_size
                Current{:,end+incmnt}={};
                cur_size=cur_size+incmnt;
            end
            Current(:,iter)={[last_ind,wrt(irt)];Next_Pointers(irt)};
        end
        Current=Current(:,1:iter);
    end
    disp([' order =',int2str(io)])
    disp(cell2mat(Current(1,:)'))
    Prev=Current;
end