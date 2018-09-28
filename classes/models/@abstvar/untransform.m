function p=untransform(self,p)
% using the linres, untransform params before evaluating the
% likelihood and later on the prior as well.
linres=self.estim_.linres;

if ~isempty(linres)
    
    p=linres.a2tilde_to_a(p);
    
end

end