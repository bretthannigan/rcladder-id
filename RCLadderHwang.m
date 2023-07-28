function [sys_id, T] = RCLadderHwang(sys)
%RCLADDERHWANG RC Ladder Reconstruction Using Routh Array Method
%   Computes the state transformation to transform an RC-ladder type system
%   to the form from RCLADDERN using the Routh array/Markov coefficient
%   algorithm from [1].
%
%   Input:
%       sys: state-space system for RC ladder system.
%
%   Outputs:
%       sys_id: strucutred system in the form of RCLadderN(R, C, 'ascending', false).
%       T_id: similarity transformation such that 
%             sys_id.A = T_id\sys.A*T_id
%
%   See also: SS2CAUERI.M, RCLADDERDIAGSCALING.M
%       [1] Hwang, C., Nakano, B., & Asahi, T. (1984). Transformation of 
%           state-space model to Cauer I and II CFE canonical forms. 
%           International Journal of Systems Science, 15(7), 797–804. 
%           https://doi.org/10.1080/00207728408926600
%
%   $Author: BH$    $Date: 2023-06-21$  $Revision: 0$ 

    n = size(sys.B, 1);
    [sys_hwang, T_hwang] = ss2caueri(sys);
    T_flip = fliplr(diag([-sys_hwang.D/sys_hwang.C(1); ones(n-1, 1)]));
    sys_id_unscaled = ss(T_flip\sys_hwang.A*T_flip, T_flip\sys_hwang.B, sys_hwang.C*T_flip, sys_hwang.D);
    [A_scaled, T_id] = RCLadderDiagScaling(sys_id_unscaled.A, sys_id_unscaled.B);
    sys_id = ss(A_scaled, sys_id_unscaled.B, sys_id_unscaled.C, sys_id_unscaled.D);
    T = T_hwang\T_flip*T_id;
end