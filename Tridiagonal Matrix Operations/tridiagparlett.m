%TRIDIAGPARLETT Trdiagonalization Procedure from Parlett1992
%   Following the algorithm and notation from [1, p. 572]. Given arbitrary 
%   square matrix B and starting vectors p1 and q1, computes similarity 
%   transformation matrices P and Q, diagonal matrix Omega, and symmetric
%   tridiagonal matrix T_hat, such that:
%       B*Q = Q/Omega*T_hat     [1, Eq. 2.1]
%       P'*B = T_hat/Omega*P'   [1, Eq. 2.2]
%   Or, alternatively [1, Theorem 2.2]:
%       P'*Q = Omega
%       P'*B*Q = Omega\T_hat
%
%   Input:
%       A: (n x n) arbitrary square matrix.
%       q1: (optional, default [1 0 ... 0]) first column vector of
%           similarity transform Q
%       p1: (optional, default [1 0 ... 0])
%           first column vector of similarity transformation P, where 
%           P = inv(Q')
%
%   Outputs:
%       T_hat: (n x n) symmetric tridiagonal matrix.
%       Q: (n x n) similarity transformation matrix.
%       Pt: (n x n) P'=Omega/Q.
%       Omega: (n x n) diagonal matrix with elements equal to the
%              codiagonal elements of T_hat.
%           
%   See also:
%       [1] Parlett, B. N. (1992). Reduction to Tridiagonal Form and 
%           Minimal Realizations. SIAM Journal on Matrix Analysis and 
%           Applications, 13(2), 567�593.
%
%   $Author: BH$    $Date: 2023-05-18$  $Revision: 0$
%
%   �2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

function [T_hat, Q, Pt, Omega] = tridiagparlett(B, q1, p1)
        
    Q = zeros(size(B));
    P = zeros(size(B));
    n = length(B);
    omega = zeros(n, 1);
    alpha = zeros(n, 1);

    P(:,1) = p1;
    Q(:,1) = q1;    

    omega(1) = P(:,1)'*Q(:,1);
    alpha(1) = P(:,1)'*B*Q(:,1);

    if length(B)>1
        Q(:,2) = B*Q(:,1) - Q(:,1)*(alpha(1)/omega(1));
        P(:,2) = P(:,1)'*B - (alpha(1)/omega(1))*P(:,1)';
        omega(2) = P(:,2)'*Q(:,2);
        alpha(2) = P(:,2)'*B*Q(:,2);
    end
    if length(B)>2
        for ii=3:length(B)
            Q(:,ii) = B*Q(:,ii-1) - Q(:,ii-1)*(alpha(ii-1)/omega(ii-1)) - Q(:,ii-2)*(omega(ii-1)/omega(ii-2));
            P(:,ii) = P(:,ii-1)'*B - (alpha(ii-1)/omega(ii-1))*P(:,ii-1)' - (omega(ii-1)/omega(ii-2))*P(:,ii-2)';
            omega(ii) = P(:,ii)'*Q(:,ii);
            if abs(omega(ii))<10*eps
                warning('The inner product p^T*q must not equal zero.')
            end
            alpha(ii) = P(:,ii)'*B*Q(:,ii);
        end
    end
    Pt = P';
    T_hat = Pt*B*Q;
    Omega = diag(omega);
    
end