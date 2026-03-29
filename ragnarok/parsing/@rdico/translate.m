%--- help for dsge/translate ---
%
%  translate -- Translates RISE codes into comprehensible atoms
% 
%  ::
% 
%    outList=translate(obj,inList)
%    outList=translate(obj,inList,order)
% 
%  Args:
% 
%     - obj (rise | dsge): scalar model object.
% 
%     - inList (char | cellstr): List of atoms to translate e.g. v_3_1, v(3,1),
%       c(4), c_4
% 
%     - order (integer|0|{1}): If order>1 the returned list of atoms is a list
%       of kroneckers in which the number of elements in each items is the
%       order of the kronecker. This is useful for instance to understand
%       what combination of variables makes up a column in a
%       differentiation. If order=0, the translation is with respect to the
%       static, rather than the dynamic model
% 
%  Returns:
% 
%     - **outList** [cellstr]: List of translated atoms
% 
%  Examples::
% 
%     translate(m,{'v(1,1)','v_1_1','c(4)','c_4'})
% 
%     list=obj.routines.symbolic_probs_times_dynamic{2}
%     outList=translate(obj,list)
%     outList=translate(obj,list,3)
%
%    Other uses of translate
%
%       polyshape/translate    rdico/translate
%