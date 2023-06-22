%TRIDIAGLANCZOS Biorthogonal Unsymmetric Lanczos Tridiagonalization
%   Following the algorithm and notation from Section 9.4.3 of [1].
%
%   Input:
%       A: (n x n) arbitrary square matrix.
%       q1: (optional, default [1 0 ... 0]) first column vector of
%           similarity transform Q
%       p1: (optional, default [1 0 ... 0])
%           first column vector of similarity transformation inverse Q^-T.
%
%   Outputs:
%       T: (n x n) nonsymmetric tridiagonal matrix.
%       Q: (n x n) similarity transformation matrix, such that T = Q\A*Q.
%       Qinv: (n x n) inverse of Q (provided for numerical stability).
%           
%   See also:
%       [1] Golub, Gene H., & Van Loan, Charles F. (1996). QRT: Matrix 
%           Computations (3rd Edition). John Hopkins University Press.
%
%   $Author: BH$    $Date: 2023-05-08$  $Revision: 1$
%
%   REVISION 1: 2023-05-15
%       Added input parser to allow specified beta vector.

function [T, Q, Qinv] = tridiaglanczos(A, varargin)
    p = inputParser;
    addRequired(p, 'A', @ismatrix);
    addOptional(p, 'p1', [], @(x) isempty(x) || isvector(x));
    addOptional(p, 'q1', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'beta', [], @isvector);
    parse(p, A, varargin{:});
    A = p.Results.A;
    p1 = p.Results.p1;
    if isempty(p1)
        p1 = [1; zeros(length(A)-1, 1)];
    else % Check normality of vectors.
        if abs(norm(p1)-1)>1e-3
            p1 = p1./norm(p1);
        end
    end
    if isrow(p1)
        p1 = p1';
    end
    q1 = p.Results.q1;
    if isempty(q1)
        q1 = [1; zeros(length(A)-1, 1)];
    else % Check normality of vectors.
        if abs(norm(q1)-1)>1e-3
            q1 = q1./norm(q1);
        end
    end
    if isrow(q1)
        q1 = q1';
    end
    if (abs(p1'*q1)<1e-6)
        error('The inner product p1''*q1 must not equal zero.')
    end
    beta_in = p.Results.beta;
    
    n = length(A);
    k = 1;
    Q = zeros(n, n+1);
    P = zeros(n, n+1);
    %Q(:,2) = q1;
    %P(:,2) = p1;
    r = q1;
    s = p1;
    gamma = zeros(1, n);
    alpha = zeros(1, n);
    beta = zeros(1, n);
    
    while ~(norm(r)<1e-6) && ~(norm(s)<1e-6) && ~(abs(s'*r)<1e-6)
        if k>=(n+1)
            break
        end
        if isempty(beta_in)
            beta(k) = norm(r);
        elseif k==1
            beta(k) = 1;
        else
            beta(k) = beta_in(k-1);
        end
        gamma(k) = s'*r/beta(k);
        Q(:,k+1) = r/beta(k);
        P(:,k+1) = s/gamma(k);
        k = k + 1;
        alpha(k) = P(:,k)'*A*Q(:,k);
        r = (A - alpha(k)*eye(n))*Q(:,k) - gamma(k-1)*Q(:,k-1);
        s = (A - alpha(k)*eye(n))'*P(:,k) - beta(k-1)*P(:,k-1);
    end
    
    if k~=(n+1)
        error(['Breakdown occurred at iteration: ' num2str(k-1)])
    end
    
    T = diag(gamma(2:end), 1) + diag(alpha(2:end)) + diag(beta(2:end), -1);
    Q = Q(:,2:end);
    P = P(:,2:end);
    Qinv = P';
end