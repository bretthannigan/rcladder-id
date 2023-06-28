function [A_scaled, T] = RCLadderDiagonalScaling(A_unscaled, B)
%RCLADDERDIAGSCALING Find Diagonal Scaling to Reconstruct Tridiagonal A
% Matrix.
%   Provided tridiagonal matrix A_unscaled with ones on the k=1 diagonal:
%                [a(1,1)    1       0     0   ..      0     ]
%                [a(2,1)  a(2,2)    1     0   ..      0     ]
%   A_unscaled = [  0     a(3,2)  a(3,3)  1   ..      0     ]
%                [  :       :       :     :   `.      :     ]
%                [  0       0       0     0   ..  a(end,end)]
%   Calculates the diagonal similarity transformation T by solving a small
%   least-squares problem:
%       [prod(t)          0                  0          ..    0   ]
%       [   0     prod(t(end:-1:2))          0          ..    0   ]
%   T = [   0             0          prod(t(end:-1:3))  ..    0   ]
%       [   :             :                  :          `.    :   ]
%       [   0             0                  0          ..  t(end)]
%   that scales the codiagonal elements to obtain [A_scaled B], where all 
%   the rows sum to zero.
%              [   a(1,1)        t(1)       0     ..      0     ]
%              [a(2,1)/t(1)     a(2,2)     t(2)   ..      0     ]
%   A_scaled = [     0       a(3,2)/t(2)  a(3,3)  ..      0     ]
%              [     :            :         :     `.      :     ]
%              [     0            0         0     ..  a(end,end)]
%
%   Inputs:
%       A_unscaled: tridiagonal matrix with ones on the upper codiagonal.
%       B: state-space B vector with only one nonzero element at the last
%       index.
%
%   Outputs:
%       A_scaled: scaled tridiagonal state-space A matrix.
%       T: similarity transformation such that A_scaled = T*A_unscaled/T
%
%   See also: RECONSTRUCTIONPARLETT.M
%
%   $Author: BH$    $Date: 2023-05-01$  $Revision: 0$

    
    b = -A_unscaled(:,end) - B;
    A = A_unscaled(:,1:end-1);
%     scaling = norm(b)/norm(A);
%     b = b./norm(b);
%     A = A./norm(A);
    scaling_factors = A\b;
    T = diag([scaling_factors; 1]); 
    A_scaled = T\A_unscaled*T;
end