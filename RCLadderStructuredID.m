function [sys_id, T] = RCLadderStructuredID(sys, Cn)
%RCLADDERSTRUCTUREDID Strucutred Identication of RC Ladder System
%
%   Inputs:
%       sys: state-space system for RC ladder system.
%       Cn: capacitance closest to readout ports.
%
%   Outputs:
%       sys_id: strucutred system in the form of RCLadderN(R, C,
%       'ascending', false).
%       T_id: similarity transformation such that 
%             sys_id.A = T_id\sys.A*T_id
%
%   See also: RCLADDERN.M, RCLADDERDIAGONALSCALING.M
%
%   $Author: BH$    $Date: 2023-06-13$  $Revision: 0$

    %% Tridiagonalize A Matrix
    n = size(sys.A, 1);
    [A_tri, Q, Pt, Omega] = tridiagparlett(sys.A, (1/sys.D*Cn).*sys.B, ((-1/sys.D).*sys.C)');
    %[A_tri, Q, Qinv, Omega] = tridiagparlett(sys.A, sys.B, sys.C');
    Tflip = fliplr(eye(n));
    A_rec = Tflip*(Omega\A_tri)*Tflip;

    %% Restore Diagonal Scaling
    B_out = Tflip*Pt*sys.B;
    %C_out = sys.C*Q*Tflip;
    C_out = sys.C*Q/Omega*Tflip;
    [A_out, T_scaling] = RCLadderDiagonalScaling(A_rec, B_out);
    sys_id = ss(A_out, B_out, C_out, sys.D);
    T = Q*Tflip*T_scaling;

%     function [A_scaled, T] = diagonal_scaling(A_unscaled, B)
%     %DIAGONAL_SCALING Find Diagonal Scaling to Reconstruct Tridiagonal A
%     % Matrix.
%     %   Provided tridiagonal matrix A_unscaled with ones on the k=1 diagonal:
%     %                [a(1,1)    1       0     0   ..      0     ]
%     %                [a(2,1)  a(2,2)    1     0   ..      0     ]
%     %   A_unscaled = [  0     a(3,2)  a(3,3)  1   ..      0     ]
%     %                [  :       :       :     :   `.      :     ]
%     %                [  0       0       0     0   ..  a(end,end)]
%     %   Calculates the diagonal similarity transformation T by solving a small
%     %   least-squares problem:
%     %       [prod(t)          0                  0          ..    0   ]
%     %       [   0     prod(t(end:-1:2))          0          ..    0   ]
%     %   T = [   0             0          prod(t(end:-1:3))  ..    0   ]
%     %       [   :             :                  :          `.    :   ]
%     %       [   0             0                  0          ..  t(end)]
%     %   that scales the codiagonal elements to obtain [A_scaled B], where all 
%     %   the rows sum to zero.
%     %              [   a(1,1)        t(1)       0     ..      0     ]
%     %              [a(2,1)/t(1)     a(2,2)     t(2)   ..      0     ]
%     %   A_scaled = [     0       a(3,2)/t(2)  a(3,3)  ..      0     ]
%     %              [     :            :         :     `.      :     ]
%     %              [     0            0         0     ..  a(end,end)]
%     %
%     %   Inputs:
%     %       A_unscaled: tridiagonal matrix with ones on the upper codiagonal.
%     %       B: state-space B vector with only one nonzero element at the last
%     %       index.
%     %
%     %   Outputs:
%     %       A_scaled: scaled tridiagonal state-space A matrix.
%     %       T: similarity transformation such that A_scaled = T*A_unscaled/T
%     %
%     %   See also: RECONSTRUCTIONPARLETT.M
%     %
%     %   $Author: BH$    $Date: 2023-05-01$  $Revision: 0$
% 
%         b = -A_unscaled(:,end) - B;
%         A = A_unscaled(:,1:end-1);
%         scaling_factors = A\b;
%         T = diag([scaling_factors; 1]); 
%         A_scaled = T\A_unscaled*T;
%     end
end