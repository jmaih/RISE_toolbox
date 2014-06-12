function [Q,the_ranks,tags]=sort_Q(Q)

the_ranks=cell2mat(Q(2,:));
the_ranks=the_ranks(:);
[the_ranks,tags]=sort(the_ranks,1,'descend');
Q=Q(:,tags);
