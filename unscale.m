function unscaledsys = unscale(sys, info)
%UNSCALE Undo the MATLAB PRESCALE command.
%   Reverse the prescaling of a state-space lti model according to the
%   equations provided in [1].
%
%   Inputs:
%       sys: scaled system.
%       info: info structure from PRESCALE command.
%
%   Outputs:
%       unscaledsys: original system.
%
%   See also: PRESCALE
%             [1] https://www.mathworks.com/help/control/ref/ss.prescale.html
%
%   $Author: BH$    $Date: 2023-06-23$  $Revision: 0$
%
%   ©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

    if ~sys.scaled
        warning('provided sys is not prescaled, continuing anyways.')
    end
    TL = diag(info.SL);
    TR = diag(info.SR);
    A = TL*sys.A*TR;
    B = TL*sys.B;
    C = sys.C*TR;
    D = sys.D;
    unscaledsys = ss(A, B, C, D, sys.Ts);
end