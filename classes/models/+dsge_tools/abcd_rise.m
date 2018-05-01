function [test,A,B,C,D]=abcd_rise(m)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% function ABCD_test(M_,options_,oo_)
%   computes ABCD Test statistics of Fernandez-Villaverde, Rubio-Ramirez,Sargent, 
%   and Watson (2007), "ABCs (and Ds) of Understanding VARs", American
%   Economic Review, 97(3), 1021-1026
%
% Algorithm
%   Uses the state-space system
%           x(+1)=Ax+Bw(+1)
%           y(+1)=Cx+Dw(+1)
%   and computes the eigenvalues of 
%           A-BC^(-1)D  (Condition 1)
%   If all eigenvalues are smaller than 1, the poor man's invertibility
%   condition is satisfied and the structural shocks can be recovered from
%   the observables
% 

if m.markov_chains.regimes_number>1
    
    error('ABCD test not developed for regime switching models')
    
end

nvobs=m.observables.number(1);

if nvobs==0
    
    error('observable variables have to be declared')
    
end

if nvobs~=sum(m.exogenous.number)
    
    error(['For this test to be performed, the number of observables has ',...
        'to be equal to the number of shocks'])
    
end

if any(m.parameters.is_measurement_error)
    
    error(['For this test to be performed, there should be no measurement ',...
        'errors'])
    
end

[m,retcode]=solve(m);

if retcode
    
    error(['model could not solve:: ',decipher(retcode)])
    
end

[T,R]=load_solution(m,'iov');

pred=m.endogenous.is_predetermined;

obs=m.observables.state_id;

A=T{1}(pred,pred);

B=R{1}(pred,:);

C=T{1}(obs,pred);

D=R{1}(obs,:);

eigvals=eig(A-B/D*C);

test=all(abs(eigvals)<1);
