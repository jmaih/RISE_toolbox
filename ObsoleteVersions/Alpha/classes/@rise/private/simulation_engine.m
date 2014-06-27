function x1=simulation_engine(solution,x0,e,sig,st,solve_order,irf_anticipate)
if nargin<6
    solve_order=[];
    if nargin<5
        st=[];
    end
end
if isempty(solve_order),solve_order=1;end

ss=solution.ss;
if isempty(st)
    if numel(ss)==1
        error('multiple regimes detected: state must be specified')
    end
    st=1;
end

% simulate zeroth order
%---------------------
x1=ss{st};
if solve_order>0
    if size(e,2)>1 && solve_order>1
        error('simulation with future information not yet implemented for higher-order approximations')
    end
    % solution matrices for order 1
    %------------------------------
    m_x=solution.m_x;
    m_e=solution.m_e;
    m_sig=solution.m_sig;
    % simulate first order
    %---------------------
    x0_hat=x0-ss{st};
    x1=x1+m_x{st}*x0_hat;
    t=0;
    while t<size(e,2)
        t=t+1;
        if irf_anticipate||t==1
            x1=x1+m_e{st}(:,:,t)*e(:,t);
        end
    end
    x1=x1+m_sig{st}*sig;
    if solve_order>1
        % solution matrices for order 2
        %------------------------------
        m_x_x=solution.m_x_x;
        m_x_e=solution.m_x_e;
        m_x_sig=solution.m_x_sig;
        m_e_e=solution.m_e_e;
        m_e_sig=solution.m_e_sig;
        m_sig_sig=solution.m_sig_sig;
        % simulate second order
        %----------------------
        x1=x1+0.5*(m_x_x{st}*kron(x0_hat,x0_hat)+...
            m_e_e{st}*kron(e(:,1),e(:,1))+...
            m_sig_sig{st}*sig^2)+...
            m_x_e{st}*kron(x0_hat,e(:,1))+...
            m_x_sig{st}*x0_hat*sig+...
            m_e_sig{st}*e(:,1)*sig;
        if solve_order>2
            error('simulation implemented only up to 2nd order approximation')
        end
    end
end