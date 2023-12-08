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
%   See also: RCLADDERN.M, RCLADDERDIAGSCALING.M
%
%   $Author: BH$    $Date: 2023-06-13$  $Revision: 0$
%
%   ©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

    %% Tridiagonalize A Matrix
    n = size(sys.A, 1);
    scaling_factor = 0.6e6; % Empirical scaling factor that works well for R's scaled to 10k-100k and C's scaled around 10p-100p.
    [A_tri, Q, Pt, Omega] = tridiagparlett((1/scaling_factor)*sys.A, (1/scaling_factor)*(1/sys.D*Cn).*sys.B, (1/scaling_factor)*((-1/sys.D).*sys.C)');
    Tflip = fliplr(eye(n));
    invOmega = diag(1./diag(Omega)); % Fix because inverse of diagonal matrix Omega is numerically unstable.
    A_rec = Tflip*invOmega*A_tri*Tflip;

    %% Restore Diagonal Scaling
    B_out = Tflip*Pt*sys.B;
    C_out = sys.C*Q*invOmega*Tflip;
    [A_out, T_scaling] = RCLadderDiagScaling(A_rec, B_out);
    sys_id = ss(A_out*scaling_factor, B_out*scaling_factor, (1/scaling_factor)*C_out, sys.D);
    T = Q*Tflip*T_scaling;
end