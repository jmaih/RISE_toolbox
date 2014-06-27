function Island=recenter(Island,lb,ub,flag)
if nargin<4
    flag=1;
end
pop=size(Island,2);
lb=lb(:,ones(pop,1));
ub=ub(:,ones(pop,1));
negs=Island<lb;
pos=Island>ub;
switch flag
    case 1
        Island(negs)=lb(negs);
        Island(pos)=ub(pos);
    case 2
        Island(negs)=min(ub(negs),2*lb(negs)-Island(negs));
        Island(pos)=max(lb(pos),2*ub(pos)-Island(pos));
    case 3
        if any(negs)
            Island(negs)=lb(negs)+(ub(negs)-lb(negs)).*rand(sum(sum(negs)),1);
        end
        if any(pos)
            Island(pos)=lb(pos)+(ub(pos)-lb(pos)).*rand(sum(sum(pos)),1);
        end
    otherwise
end
end

