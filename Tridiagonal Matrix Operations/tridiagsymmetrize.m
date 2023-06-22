%TRIDIAGSYMMETRIZE Compute the similarity transformation to symmetrize a 
%tridiagonal matrix.
%   Following along from Slide 3-4 of "Tridiagonal Matrices" by Gérard 
%   Meurant (2008) [1]. Assumes that the upper codiagonal coefficients 
%   omega_j are nonzero and that all products omega_j*beta_j are positive, 
%   where omega is the lower codiagonal.
%
%   Input:
%       T: (n x n) tridiagonal matrix satisfying the conditions above.
%
%   Outputs:
%       Tsym: (n x n) symmetric tridiagonal matrix similar to T.
%       D: (n x n) similarity transformation matrix such that Tsym = D*T/D.
%           
%   See also:
%       [1] Meurant, Gerard. Tridiagonal Matrices (2008). 
%           URL: https://www.math.hkbu.edu.hk/ICM/LecturesAndSeminars/08OctMaterials/1/Slide3.pdf
%           [Accessed: 2023-05-02]
%       [2] Maple script "TridiagonalMatrixSymmetrization.mw"
%
%   $Author: BH$    $Date: 2023-05-01$  $Revision: 0$

function [Tsym, D] = tridiagsymmetrize(T)

    n = length(T);
    omega = diag(T, 1);
    beta = diag(T, -1);
    if any(abs(omega)<100*eps)
        error('All upper codiagonal elements of T, diag(T, 1), must be nonzero.')
    end
    if any(beta.*omega<0)
        error('All element-wise products of the upper and lower codiagonals, diag(T, 1).*diag(T, -1), must be positive.')
    end

    delta = zeros(n, 1);
    delta(1) = 1;
    for i_delta=2:n
        delta(i_delta) = sqrt(prod(beta(1:i_delta-1))/prod(omega(1:i_delta-1)));
    end
    D = diag(delta);
    Tsym = D\T*D;
        
end