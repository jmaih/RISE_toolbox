function varargout=output_in_original_zha_matlab_code(m,scale)

nout=nargout;

if nout==0
    disp('----------------------------------------------------------------------------------------------------------------');
    disp('           Printing out A0hat, A1hat, and A2hat in a form compatible with Zha''s original Matlab code           ');
    disp('----------------------------------------------------------------------------------------------------------------');
end
Ahat=cell(1,m.nlags+1);
for ireg=1:m.markov_chains.regimes_number
    if nout==0
        fprintf('********** Regime # %0.0f ***********\n',ireg);
    end
    for ilag=0:m.nlags
        if nout==0
            fprintf('A%0.0fhat \n',ilag);
        end
        lagname=sprintf('a%0.0f',ilag);
        tmp = ((m.solution.sig{ireg}\m.solution.(lagname){ireg})*scale)';
        if nout==0
            disp(tmp)
        end
        Ahat{ilag+1}(:,:,ireg)=tmp;
    end
end
if nargout
    varargout=Ahat(1:nargout);
end
