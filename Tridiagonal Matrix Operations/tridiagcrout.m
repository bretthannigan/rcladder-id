function [L, U] = tridiagcrout(A)
%CROUT Crout Factorization of a Tridiagonal Matrix
%
%   Input:
%       A: (n x n) tridiagonal matrix.
%   Outputs:
%       L: (n x n) lower trianguar matrix.
%       U: (n x n) unit upper triangular matrix.
%
%   See also: http://www2.lawrence.edu/fast/GREGGJ/MATH420/Chapt6/Section_6_6.pdf
%             https://www.webpages.uidaho.edu/~barannyk/Teaching/LU_factorization_tridiagonal.pdf
%
%   ©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

    if ~(all(all(abs(triu(A, 2))<0.01*norm(A))) && all(all(abs(tril(A, -2))<0.01*norm(A))))
        warning('Matrix A must be tridiagonal. Assuming zero entries outside the -1, 0, 1 diagonals.');
    end

    Ldiag = diag(A, -1);
    L = diag(Ldiag, -1); % l(i, i-1) = a(i, i-1)
    Lii = zeros(size(L, 1), 1);
    Uii = zeros(size(L, 1), 1);
    Lii(1) = A(1,1);
    for ii=2:length(Lii)
        Uii(ii) = A(ii-1, ii)/Lii(ii-1);
        Lii(ii) = A(ii, ii) - A(ii, ii-1)*Uii(ii);
    end
    L = L + diag(Lii);
    U = eye(size(A)) + diag(Uii(2:end), 1);

end