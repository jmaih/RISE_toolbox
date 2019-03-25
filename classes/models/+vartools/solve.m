% SOLVE  Symbolic solution of algebraic equations.
%    S = SOLVE(eqn1,eqn2,...,eqnM,var1,var2,...,varN)
%    S = SOLVE(eqn1,eqn2,...,eqnM,var1,var2,...,varN,'ReturnConditions',true)
% 
%    [S1,...,SN] = SOLVE(eqn1,eqn2,...,eqnM,var1,var2,...,varN)
%    [S1,...,SN,params,conds] = SOLVE(eqn1,...,eqnM,var1,var2,...,varN,'ReturnConditions',true)
% 
%    The eqns are symbolic expressions, equations, or inequalities.  The
%    vars are symbolic variables specifying the unknown variables.
%    If the expressions are not equations or inequalities, 
%    SOLVE seeks zeros of the expressions.
%    Otherwise SOLVE seeks solutions.
%    If not specified, the unknowns in the system are determined by SYMVAR,
%    such that their number equals the number of equations.
%    If no analytical solution is found, a numeric solution is attempted;
%    in this case, a warning is printed. 
% 
%    Three different types of output are possible.  For one variable and one
%    output, the resulting solution is returned, with multiple solutions to
%    a nonlinear equation in a symbolic vector.  For several variables and
%    several outputs, the results are sorted in the same order as the
%    variables var1,var2,...,varN in the call to SOLVE.  In case no variables
%    are given in the call to SOLVE, the results are sorted in lexicographic
%    order and assigned to the outputs.  For several variables and a single
%    output, a structure containing the solutions is returned.
% 
%    SOLVE(...,'ReturnConditions', VAL) controls whether SOLVE should in  
%    addition return a vector of all newly generated parameters to express 
%    infinite solution sets and about conditions on the input parameters 
%    under which the solutions are correct. 
%    If VAL is TRUE, parameters and conditions are assigned to the last two 
%    outputs. Thus, if you provide several outputs, their number must equal 
%    the number of specified variables plus two.
%    If you provide a single output, a structure is returned 
%    that contains two additional fields 'parameters' and 'conditions'.
%    No numeric solution is attempted even if no analytical solution is found.
%    If VAL is FALSE, then SOLVE may warn about newly generated parameters or 
%    replace them automatically by admissible values. It may also fall back 
%    to the numerical solver.
%    The default is FALSE.
% 
%    SOLVE(...,'IgnoreAnalyticConstraints',VAL) controls the level of
%    mathematical rigor to use on the analytical constraints of the solution
%    (branch cuts, division by zero, etc). The options for VAL are TRUE or
%    FALSE. Specify FALSE to use the highest level of mathematical rigor
%    in finding any solutions. The default is FALSE.
% 
%    SOLVE(...,'PrincipalValue',VAL) controls whether SOLVE should return multiple
%    solutions (if VAL is FALSE), or just a single solution (when VAL is TRUE).
%    The default is FALSE.
% 
%    SOLVE(...,'IgnoreProperties',VAL) controls if SOLVE should take
%    assumptions on variables into account. VAL can be TRUE or FALSE.
%    The default is FALSE (i.e., take assumptions into account).
% 
%    SOLVE(...,'Real',VAL) allows to put the solver into "real mode."
%    In "real mode," only real solutions such that all intermediate values
%    of the input expression are real are searched. VAL can be TRUE or FALSE.
%    The default is FALSE.
% 
%    SOLVE(...,'MaxDegree',n) controls the maximum degree of polynomials
%    for which explicit formulas will be used during the computation.
%    n must be a positive integer. The default is 3.
% 
%    Example 1:
%       syms p x r
%       solve(p*sin(x) == r) chooses 'x' as the unknown and returns
% 
%         ans =
%                asin(r/p)
%           pi - asin(r/p)
% 
%    Example 2:
%       syms x y
%       [Sx,Sy] = solve(x^2 + x*y + y == 3,x^2 - 4*x + 3 == 0) returns
% 
%         Sx =
%          1
%          3
% 
%         Sy =
%             1
%          -3/2
% 
%    Example 3:
%       syms x y
%       S = solve(x^2*y^2 - 2*x - 1 == 0,x^2 - y^2 - 1 == 0) returns
%       the solutions in a structure.
% 
%         S =
%           x: [8x1 sym]
%           y: [8x1 sym]
% 
%    Example 4:
%       syms a u v
%       [Su,Sv] = solve(a*u^2 + v^2 == 0,u - v == 1) regards 'a' as a
%       parameter and solves the two equations for u and v.
% 
%    Example 5:
%       syms a u v w
%       S = solve(a*u^2 + v^2,u - v == 1,a,u) regards 'v' as a
%       parameter, solves the two equations, and returns S.a and S.u.
% 
%       When assigning the result to several outputs, the order in which
%       the result is returned depends on the order in which the variables
%       are given in the call to solve:
%       [U,V] = solve(u + v,u - v == 1, u, v) assigns the value for u to U
%       and the value for v to V. In contrast to that
%       [U,V] = solve(u + v,u - v == 1, v, u) assigns the value for v to U
%       and the value of u to V.
% 
%    Example 6:
%       syms a u v
%       [Sa,Su,Sv] = solve(a*u^2 + v^2,u - v == 1,a^2 - 5*a + 6) solves
%       the three equations for a, u and v.
% 
%    Example 7:
%       syms x
%       S = solve(x^(5/2) == 8^(sym(10/3))) returns all three complex solutions:
% 
%         S =
%                                                         16
%          - 4*5^(1/2) - 4 + 4*2^(1/2)*(5 - 5^(1/2))^(1/2)*i
%          - 4*5^(1/2) - 4 - 4*2^(1/2)*(5 - 5^(1/2))^(1/2)*i
% 
%    Example 8:
%       syms x
%       S = solve(x^(5/2) == 8^(sym(10/3)), 'PrincipalValue', true)
%       selects one of these:
% 
%         S =
%         - 4*5^(1/2) - 4 + 4*2^(1/2)*(5 - 5^(1/2))^(1/2)*i
% 
%    Example 9:
%       syms x
%       S = solve(x^(5/2) == 8^(sym(10/3)), 'IgnoreAnalyticConstraints', true)
%       ignores branch cuts during internal simplifications and, in this case,
%       also returns only one solution:
% 
%         S =
%         16
% 
%    Example 10:
%       syms x
%       S = solve(sin(x) == 0) returns 0
%       
%       S = solve(sin(x) == 0, 'ReturnConditions', true) returns a structure expressing
%       the full solution:
% 
%       S.x = k*pi
%       S.parameters = k
%       S.conditions = in(k, 'integer')
% 
%    Example 11:
%       syms x y real 
%       [S, params, conditions] = solve(x^(1/2) = y, x, 'ReturnConditions', true)
%       assigns solution, parameters and conditions to the outputs. 
%       In this example, no new parameters are needed to express the solution:
%  
%       S = 
%       y^2
% 
%       params =
%       Empty sym: 1-by-0
%    
%       conditions =
%       0 <= y
% 
%    Example 12:
%       syms a x y
%       [x0, y0, params, conditions] = solve(x^2+y, x, y, 'ReturnConditions', true)
%       generates a new parameter z to express the infinitely many solutions.
%       This z can be any complex number, both solutions are valid without 
%       restricting conditions:
%       
%       x0 =
%       -(-z)^(1/2)
%       (-z)^(1/2)
% 
%       y0 =
%       z
%       z
% 
%       params =
%       z
% 
%       conditions =
%       true
%       true
% 
%    Example 13:
%       syms t positive
%       solve(t^2-1)
% 
%         ans =
%         1
% 
%       solve(t^2-1, 'IgnoreProperties', true)
% 
%         ans =
%           1
%          -1
% 
%    Example 14:
%       solve(x^3-1) returns all three complex roots:
% 
%         ans =
%                              1
%          - 1/2 + (3^(1/2)*i)/2
%          - 1/2 - (3^(1/2)*i)/2
% 
%       solve(x^3-1, 'Real', true) only returns the real root:
% 
%         ans =
%         1
% 
%    See also DSOLVE, SUBS.
%
%    Reference page in Doc Center
%       doc solve
%
%    Other functions named solve
%
%       abstvar/solve    dsge/solve    rfvar/solve
%