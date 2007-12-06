%  species tree
%       4
%      / \
%     3   \
%    / \   \
%   /   \   \
%   0   1    2
%   A   B    C    # species names
%
%  gene tree
%
%
%           10
%          /\ 
%         9  \
%        /\   \
%       /  \   \
%      8    \   \ 
%     / \    \   \
%    6   \    7   \
%   / \   \  / \   \
%  0  1   2 3  4   5
%  A  A   B A  B   C
% 

function test_reconcile()
    nnodes = 11
    nsnodes = 5

    % gene parent tree (ptree) length=2n-1
    ptree = [6, 6, 8, 7, 7, 10, 8, 9, 9, 10, -1]
    
    % species parent tree (pstree)  length=2m-1
    pstree = [3, 3, 4, 4, -1]
    
    % map gene leaves to species leaves (gene2species) length=2n-1
    % gene2species[gene node's id] = species node's id
    gene2species = [0, 0, 1, 0, 1, 2, -1, -1, -1, -1, -1]
    
    
    
    % display gene tree for debuging
    spidir_display_tree(ptree, ones(1, nnodes), 10)
    spidir_display_tree(pstree, ones(1, nsnodes), 10)
    
    % make the call
    [recon, events] = spidir_reconcile(ptree, pstree, gene2species)
    
    % for each gene node i
    % events[i] = 0   if i is a gene
    %           = 1   if i is a speciation
    %           = 2   if i is a duplication
    
    % for each gene node i
    % recon[i] = j
    % where j is a species node
