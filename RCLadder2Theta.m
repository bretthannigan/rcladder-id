function [R, C] = RCLadder2Theta(sys)
%RCLADDER2THETA Extract R and C parameters from tridiagonal state-space
% form.
%
%   Inputs:
%       sys: state-space system for RC ladder system, already placed in the 
%            correct form of RCLADDERN(R, C, 'ascending', false)
%
%   Outputs:
%       R: vector of resistances.
%       C: vector of capacitances.
%
%   See also: RCLADDERN.M, TRIDIAGCROUT.M
%
%   $Author: BH$    $Date: 2023-06-21$  $Revision: 0$

    ABCD = [sys.A, sys.B; sys.C, sys.D];
    [U, L] = tridiagcrout(ABCD');
    C = -cumprod(flipud(diag(L', -1))).*cumprod(-ones(length(sys.A), 1));
    R = 1./(diag(U, -1).*flipud(C));
end