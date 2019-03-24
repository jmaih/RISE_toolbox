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
%     obj (rise | dsge): scalar model object.
% 
%     inList (char | cellstr): List of atoms to translate e.g. y_3,
%        param(4), ss_20, def_10
% 
%     order (integer|0|{1}): If order>1 the returned list of atoms is a list
%        of kroneckers in which the number of elements in each items is the
%        order of the kronecker. This is useful for instance to understand
%        what combination of variables makes up a column in a
%        differentiation. If order=0, the translation is with respect to the
%        static, rather than the dynamic model
% 
%  Returns:
%     :
% 
%     - **outList** [cellstr]: List of translated atoms
% 
%  Example:
%     :
% 
%     list=obj.routines.symbolic.probs_times_dynamic{2}
%     outList=translate(obj,list)
%     outList=translate(obj,list,3)
% 
%
%    Other functions named translate
%
%       polyshape/translate
%