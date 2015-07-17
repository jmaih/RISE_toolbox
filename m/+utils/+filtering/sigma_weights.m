function [w]=sigma_weights(m,type)

switch type
    case 'ckf'
        % 2m points
        w = 1/(2*m)*ones(1,2*m);
    case 'ddf'
        % 2m+1 points
        delta=sqrt(3);
        w = [(delta^2-m)/(delta^2),1/(2*delta^2)*ones(1,2*m)];
    case 'ukf'
        % 2m+1 points
        k=3-m;
        w=1/(2*(m+k))*ones(1,2*m+1);
        w(1)=k/(m+k);
    otherwise
        error(['unknown type of rule ',type])
end

end
