    % Substract from the data y(t), the mean Z x B(t)
    yhat = zeros(M,t);
    for i = 1:t
        yhat(:,i) = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw(:,i);
    end
    
    % This part is more tricky, check Primiceri
    % Zc is a [M x M(M-1)/2] matrix defined in (A.2) page 845, Primiceri
    Zc = - yhat(:,:)';
    for i = 1:t
        sigmatemp = diag(sigt((i-1)*M+1:i*M,:)^2)';
        sigma2temp(i,:) = sigmatemp;
    end
    
    Atdraw = [];
    ind = 1;
    for ii = 2:M
        % Draw each block of A(t)
        [Atblockdraw,log_lik2a] = carter_kohn(yhat(ii,:),Zc(:,1:ii-1),sigma2temp(:,ii),...
            Sblockdraw{ii-1},sizeS(ii-1),1,t,A_0_prmean(((ii-1)+(ii-3)*(ii-2)/2):ind,:),A_0_prvar(((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind));
        Atdraw = [Atdraw ; Atblockdraw]; %#ok<AGROW> % Atdraw is the final matrix of draws of A(t)
        ind = ind + ii;
    end
    
    %=====| Draw S, the covariance of A(t) (from iWishart)
    % Take the SSE in the state equation of A(t)
    Attemp = Atdraw(:,2:t)' - Atdraw(:,1:t-1)';
    sse_2 = zeros(numa,numa);
    for i = 1:t-1
        sse_2 = sse_2 + Attemp(i,:)'*Attemp(i,:);
    end
    % ...and subsequently draw S, the covariance matrix of A(t) 
    ijc = 1;
    for jj=2:M
        Sinv = inv(sse_2(((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc) + S_prmean{jj-1});
        Sinvblockdraw = wish(Sinv,t+S_prvar(jj-1));
        Sblockdraw{jj-1} = inv(Sinvblockdraw); % this is a draw from S
        ijc = ijc + jj;
    end