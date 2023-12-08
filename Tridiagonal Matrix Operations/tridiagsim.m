%TRIDIAGSIM Tridiagonalization via Similarity Transformations
%   Following the "Elementary Similarity Transformations" method from
%   Section 2.1 of [1]. Recursive implementation. Probably suffers from
%   breakdown in certain cases. This is a special (simple) case of Lanczos
%   tridiagonalization [2] with starting vectors u_1 = e_1 and 
%   v_1 = e_1.
%
%   Input:
%       T: (n x n) arbitrary square matrix.
%
%   Outputs:
%       A: (n x n) nonsymmetric tridiagonal matrix.
%       S: (n x n) similarity transformation matrix, such that A = S*T/S.
%       Sinv: (n x n) inverse of S (provided for numerical stability).
%           
%   See also:
%       [1] Sidje, R. B., & Burrage, K. (2005). QRT: A QR-Based 
%           Tridiagonalization Algorithm for Nonsymmetric Matrices. SIAM 
%           Journal on Matrix Analysis and Applications, 26(3), 878–900. 
%           https://doi.org/10.1137/040612476
%       [2] TRIDIAGLANCZOS.M
%
%   $Author: BH$    $Date: 2023-04-27$  $Revision: 0$
%
%   ©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

function [A, S, Sinv] = tridiagsim(T)
    
    delta = T(1,1);
    if length(T)>1 % Recursive case.
        n = length(T);
        e1 = [1; zeros(n-2, 1)];
        y = T(1,2:end)';
        x = T(2:end,1);
        Z = T(2:end,2:end);

        M2 = eye(n-1) - 1/x(1)*x.*e1' - 1/y(1)*e1*y';
        M2inv = eye(n-1) - e1*e1' - (y'*x)\(x*y');

        [A_red, S_red, Sinv_red] = tridiagsim(M2*Z*M2inv);
        A = [delta y'*M2inv; M2*x A_red];
        % Originally thought the multiplication order in the next 2 lines
        % should be reversed.
        S = blkdiag(1, S_red)*blkdiag(1, M2);
        Sinv = blkdiag(1, M2inv)*blkdiag(1, Sinv_red);
    else % Base case.
        A = delta;
        S = 1;
        Sinv = 1;
    end
    
end