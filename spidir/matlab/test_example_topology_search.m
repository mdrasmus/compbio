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


function test_example_topology_search()
    nnodes = 7
    
    % example input sequence
    seqs = ['AAAAGGGAAAAAAA';
            'GGGGGGGAAAAAAA';
            'CCCCCCCAAAAAAA';
            'TTTTCCCAAAAAAA']
        
    % defaults for ML HKY DIST
    bgfreq = [.25, .25, .25, .25]
    tsvratio = .5    
    maxiter = 10
    
    % initial proposed topology from NJ
    distmat = align2distmat(seqs)
    [ptree, dists] = spidir_neighborjoin(distmat)
    spidir_display_tree(ptree, dists, 80.0)
        
    % start a search loop
    maxlogl = -inf;
    maxtop = ptree;
    
    for i = 1:10
        % evaluate the likelihood of the current topology
        [logl, dists] = spidir_mlhkydist(ptree, seqs, bgfreq, tsvratio, maxiter)
        spidir_display_tree(ptree, dists, 20)
        
        % is likelihood the highest we have seen?
        if logl > maxlogl
            maxlogl = logl;
            maxtop = ptree;
        end
        
        % change topology
        % this is where you do your SPR/NNI/reconciliation/etc
    end
    
    % display best tree
    maxtop
    maxlogl


% make a distance matrix from an alignment
function distmat = align2distmat(align)
    [nseqs, seqlength] = size(align)
    distmat = zeros(nseqs)
    
    for i = 1:nseqs
        for j = 1:nseqs
            changes = 0
            
            for k = 1:seqlength
                if align(i+(k-1)*nseqs) ~=  align(j+(k-1)*nseqs)
                    changes =  changes + 1
                end
            end
            
            distmat(i, j) = changes / seqlength
        end
    end
    
