function varargout=integrate_regimes(Q,varargin)
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

reg_nbr=size(Q,1);

assert(size(Q,2)==reg_nbr,'transition matrix should be square')

assert(all(abs(sum(Q,2)-1)<sqrt(eps)),'Bad transition matrix')

nout=nargout;

assert(nout==length(varargin),...
    'number of output arguments must be equal to the number of input arguments minus 1')

ergodic=true;

[pai,retcode]=initial_markov_distribution(Q,ergodic);

if retcode
    
    warning(decipher(retcode))
    
end

varargout=num2cell(zeros(1,nout));

for ireg=1:reg_nbr
    
    for iout=1:nout
        
        if ireg==1 
            
            if ~iscell(varargin{iout})
            
                error('all input arguments after the first one should be cell arrays')
            
            end
            
            if reg_nbr~=numel(varargin{iout})
                
                error('all input arguments after the first should match the number of regimes')
            
            end
            
        end
        
        varargout{iout}=varargout{iout}+pai(ireg)*varargin{iout}{ireg};
    
    end
    
end

end