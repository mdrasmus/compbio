%
%         6  
%        / \
%       /   \
%      /     \
%     4       5
%    / \     / \
%   /   \   /   \
%   0   1   2   3
%


function test_mlhkydist()
    nnodes = 4
    nsnodes = 4

    % gene parent tree (ptree) length=2n-1
    ptree = [4, 4, 5, 5, 6, 6, -1]
    
    seqs = ['AAAAGGGAAAAAAA';
            'GGGGGGGAAAAAAA';
            'CCCCCCCAAAAAAA';
            'TTTTCCCAAAAAAA']
    
    bgfreq = [.25, .25, .25, .25]

    tsvratio = .5
    
    maxiter = 10
    
    
    % make the call
    %[logl, dists] = spidir_mlhkydist(ptree, seqs, bgfreq, tsvratio, maxiter)
    %spidir_display_tree(ptree, dists, 20)
    
    nnodes = 4
    nsnodes = 4

    % gene parent tree (ptree) length=2n-1
    ptree = [4, 4, 5, 5, 6, 6, -1]

    seqs = ['AAAAGGGAAAAAAA';
            'GGGGGGGAAAAAAA';
            'CCCCCCCAAAAAAA';
            'TTTTCCCAAAAAAA']

    bgfreq = [.25, .25, .25, .25]

    tsvratio = .5

    maxiter = 10
    
    
    % make the call
    [logl, dists] = spidir_mlhkydist(ptree, seqs, bgfreq, tsvratio, maxiter)
    spidir_display_tree(ptree, dists, 100)