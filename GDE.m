function [D_rank_A] = GDE(M)
% April 2021
%
% This matlab code implements the Rank Estimation of Low-Rank Matrix for ARE-RPCA in 
% "Adaptive Rank Estimate in Robust Principal Component Analysis".
% Xu Z, He R, Xie S, et al.
%
% M - m x n matrix of observations/data (required input)
%
% D_rank_A - The estimated rank 
%    

    num = 1 ; % The adjustment parameter of D_M

    [m,n] = size(M);
    if m > n
        R_M = M' * M;
    else
        R_M = M * M';
    end
    RR = R_M;
    [V_eig, D_eig] = eig(RR);
    % Order the eigenvalues
    [D_sort,index] = sort(diag(D_eig),'descend');
    D_sort = diag(D_sort(index));
    V_sort = V_eig(:,index);
    
    % Reduction correlation matrix
    R = V_sort*D_sort*V_sort';

    [M1,N1] = size(R);
    for p1=1:M1-1
        for q1=1:M1-1
            R_A_T(p1,q1)=R(p1,q1);
        end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    end
    [V,D] = eig(R_A_T);
    t = zeros(M1-1,1);
    T = [V t;t' 1];

    R1 = R;
    T1 = T;
    R_T_A = zeros(N1,N1);
    rA = zeros(1,N1-1);
    R_T = T1'*R1*T1;
    
    % The radii of the Gerschgorin¡¯s disk
    for i=1:N1-1
        r_r_1(i)=abs(R_T(i,N1));
    end

    D1=diag(D);  

    % The shrunk paremeter
    D_seta_N=1/sqrt(D1'*D1);
    D_seta = D1*D_seta_N;
    % The shrunk radius
    r_r = r_r_1';
    r_r = r_r.*D_seta;
    r=abs(r_r');

    for i=1:N1-2
        D_M=abs(2*D1(i+1))/(num*sqrt(D1(i:N1-1)'*D1(i:N1-1))); 
        m_seta_r=sum(r);
        D_MM = D_M*(m_seta_r/(N1-1));
        % The improved heuristic decision rule
        GDE_K(i)=r(i)-D_MM; 
        if GDE_K(i)<0 
             D_rank_A=i-1;
             break;
        end
    end
end
