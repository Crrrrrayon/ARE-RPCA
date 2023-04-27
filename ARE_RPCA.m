function [L, S, iter, psnr_n] = ARE_RPCA(D, A)
%% Rank N soft constraint for RPCA (Ver. 1.0)
%
%	This code follows and refers the Copyright and Intelligent Property codes of Wuhan University of Science and Technology.
%	<Input>
%	D : [m x n] observation matrix (Notice: must be m>>n)
%	rankN : target rank N
%
%	<Output>
%	L : [m x n] low-rank matrix
%	S : [m x n] sparse matrix
%
%%
%
% This matlab code implements the Adaptive RPCA based on Iterative Rank Estimate in
% "Adaptive Rank Estimate in Robust Principal Component Analysis".
% Xu Z, He R, Xie S, et al.
%
%   This release is for research usage.
%   Email me with questions and comments: xuzhengqin@wust.edu.cn, I'll help you.
%
%                                                                   Seoul, Korea.
%                                                                   Jan.
%                                                                   3. 2014.
%
%%
[m n] = size(D);
lambda = 1 / sqrt(m);
tol = 1e-7;
maxIter = 1000;
% initialize
Y = D;						% can be tuned
norm_two = norm(Y, 2);
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
L = zeros( m, n);			% can be tuned
S = zeros( m, n);			% can be tuned
mu = 1.25/norm_two; 		% can be tuned
mu_bar = mu * 1e7;
rho = 1.05;          		% can be tuned
d_norm = norm(D, 'fro');
rank_N(1) = GDE(D);
rankN=rank_N(1);
% Initial value
k_r = 0;
n_r = 0;

iter = 0;
converged = false;
stopCriterion = 1;
index=1;
while ~converged
    iter = iter + 1;
    
    %% E
    temp_T = D - L + (1/mu)*Y;
    S = max(temp_T - lambda/mu, 0);
    s1 = temp_T + lambda/mu;
    s2 = min(temp_T + lambda/mu, 0);
    s3 = length(find(s2<0));
    S = S+min(temp_T + lambda/mu, 0);
    
    %% A
    [U, Singular, V] = svd(D - S + (1/mu)*Y, 'econ');
    diagS = diag(Singular);
    svp = length(find(diagS > 1/mu));
    
    if svp <= rankN
        dterm = diag( diagS(1:svp) );
    else
        dterm = diag( [diagS(1:rankN);diagS(rankN+1:svp) - 1/mu] );
    end
    L = U(:, 1:svp) * dterm * V(:, 1:svp)';
    
    if rem(iter,2) == 1
        psnr_n(index) = psnr(A ,L);
        disp(['GDE_PSNR(',num2str(index),')===',num2str(psnr_n(index))]);
        index = index+1;
    end
    
    
    if k_r ~ = 4
        rak(iter) = rank(L); % The rank of the low-rank matrix
        rank_N(iter+1) = GDE(L);
        rankN = rank_N(iter+1);
        if rankN == 0
            rankN = rak(iter);
        end
        n_r = n_r +1;
        if rak(iter) == rank_N(iter+1)
            k_r = k_r + 1;
        else
            k_r = 0;
        end
    else
        rankN = rank(L);
        disp(['Rank =',num2str(rankN),',iter = ',num2str(n_r)])
    end
    
    Z = D - L - S;
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
        disp(stopCriterion);
    end
    
    if  iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
    if iter<1000
        psnr_n(index:500)=psnr_n(index-1);
    end
end
end

