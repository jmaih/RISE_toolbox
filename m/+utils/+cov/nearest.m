%--- help for digraph/nearest ---
%
%  NEAREST Compute nearest neighbors of a node
% 
%    NODEIDS = NEAREST(G, S, D) returns all nodes within a distance D
%    from node S, sorted from nearest to furthest.  If the graph is
%    weighted (that is G.Edges contains a Weight variable) those
%    weights are used as the distances along the edges in the graph.
%    Otherwise, all distances are implicitly taken to be 1.
% 
%    [NODEIDS, DIST] = NEAREST(G, S, D) additionally returns a
%    vector DIST, where DIST(j) is the distance from node SOURCE to node
%    NODEIDS(j).
% 
%    NEAREST(..., 'Direction', DIR) specifies the search direction. DIR can
%    be:
%        'outgoing' - Distances are computed from node S to nodes NODEIDS
%                     (this is the default).
%        'incoming' - Distances are computed from nodes NODEIDS to node S.
% 
%    NEAREST(...,'Method',METHODFLAG) optionally specifies the method to
%    compute the distances.
%    METHODFLAG can be:
% 
%          'auto'  -  Uses 'unweighted' if no weights are set, 'positive'
%                     if all weights are nonnegative, and 'mixed' otherwise.
%                     This method is the default.
% 
%    'unweighted'  -  Treats all edge weights as 1.
% 
%      'positive'  -  Requires all edge weights to be positive.
% 
%         'mixed'  -  Allows negative edge weights, but requires that the
%                     graph has no negative cycles.
% 
%    See also SHORTESTPATHTREE, DISTANCES, GRAPH/NEAREST
%
%    Reference page in Doc Center
%       doc digraph/nearest
%
%    Other functions named nearest
%
%       graph/nearest
%