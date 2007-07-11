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

function test_matlab_spidir()
    nnodes = 4
    nsnodes = 4

    % gene parent tree (ptree)
    ptree = [3, 3, 4, 4, -1]
    
    % gene branch lengths (dists)
    dists = [1, 1, .5, 2, 0]
    
    % species parent tree (pstree)
    pstree = [3, 3, 4, 4, -1]
    
    % map gene leaves to species leaves (gene2species)
    gene2species = [0, 1, 2, -1, -1]
    
    % mean and sdev of relative branch length distributions (mu, sigma)
    mu = [.9, .8, .3, 1.8, 0]
    sigma = [.2, .2, .2, .2, 0]
    
    % gene rate distribution
    alpha = 9.50454
    beta  = 5.48942
    
    % gene rate (generate)
    % specify a negative number to request SPIDIR to estimate it for you
    generate = 1.0
    
    % make the call
    logl = matlab_spidir(ptree, dists, pstree, gene2species, ...
                         mu, sigma, alpha, beta, generate)
    
    
    % compute by hand to compare
    logl2 = 0;
    for i = 1:nnodes
        lk = normal(dists(i), mu(i), sigma(i));
        logl2 = logl2 + log(lk);
    end
    logl2 = logl2 + log(gammadist(generate, alpha, beta));
    logl2
    

function g = gammadist(x, alpha, beta)
    g = (exp(-x * beta) * (x^(alpha - 1)) * (beta^alpha)) / gamma(alpha);
    

function n = normal(x, u, s)
    n = 1/(s*sqrt(2*pi)) * exp(- (x - u)^2 / (2*s*s));
