function V=derivative_indices(number_of_variables,order)
% this function returns the indices for the derivatives calculated by
% mutivariate_taylor_approximation
V=cell(1,order);
d=0;
while d<order
    d=d+1;
    new=update_index(V{d},number_of_variables,d);
    while ~isempty(new)
        V{d}=[V{d};new];
        new=update_index(new,number_of_variables,d);
    end
end

function new=update_index(old,number_of_variables,d)
if isempty(old)
    new=ones(1,d);
else
    last=old(end);
    if last<number_of_variables
        new=old;
        new(end)=new(end)+1;
    elseif last==number_of_variables
        % find the first index that is less than number_of_variables
        notfound=true;
        ii=d;
        while notfound && ii>1
            ii=ii-1;
            if old(ii)<number_of_variables
                new=old;
                new(ii:end)=new(ii)+1;
                notfound=false;
            end
        end
        if notfound
            new=[];
        end
    end
end

% function V=derivative_indices(number_of_variables,order)
% V=cell(1,order);
% for d=1:order
% %     disp(['order=',int2str(order)])
%     new=ones(1,d);
%     last=new(end);
%     V{d}=[V{d};new];
%     while last<number_of_variables
%         last=last+1;
%         new(end)=last;
%         V{d}=[V{d};new];
%     end
%     ii=d;
%     while ~all(new==number_of_variables)
%         if new(ii)<number_of_variables
%             new(ii:end)=new(ii)+1;
%             last=new(end);
%             V{d}=[V{d};new];
%             while last<number_of_variables
%                 last=last+1;
%                 new(end)=last;
%                 V{d}=[V{d};new];
%             end
%         else
%             ii=ii-1;
%         end
%     end
% end

