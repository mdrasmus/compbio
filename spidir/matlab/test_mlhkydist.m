%
% Let's say both the gene tree and species tree have this same structure
%
%       4
%      / \
%     3   \
%    / \   \
%   /   \   \
%   0   1    2
%
% Then the parent tree array representation is
%   ptree = [3, 3, 4, 4, -1]
% such that ptree[node's id] = node's parent's id
%
% In addition, the following must be true
% 1. tree must be binary
% 2. leaves must be numbered 0 to n-1
% 3. internal nodes are numbered n to 2n-1
% 4. root must be numbered 2n-1
% 5. the parent of root is -1
%

function test_mlhkydist()
    nnodes = 4
    nsnodes = 4

    % gene parent tree (ptree) length=2n-1
    ptree = [3, 3, 4, 4, -1]
    
    seqs = ['GAAAAAAA';
            'CCCAAAAA';
            'CCCTTAAA']
    
    bgfreq = [.25, .25, .25, .25]

    tsvratio = .5
    
    maxiter = 10
    
    
    % make the call
    [logl, dists] = spidir_mlhkydist(ptree, seqs, bgfreq, tsvratio, maxiter)
    spidir_display_tree(ptree, dists, 2)
    2
