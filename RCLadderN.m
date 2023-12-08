function varargout = RCLadderN(R, C, varargin)
% RCLADDERN Generate model of the admittance state-space model of an 
%   N-section RC ladder (with R_1 in last position).
%
%    i_N(t)               i_n(t)                i_1(t)
%   o->-[  R_N  ]---o---...->-[ R_(n) ]---o---...->-[  R_1  ]--,
%   +            +  |                  +  |                  + |
%                   |                     |                    |
%  v_in(t)  v_N(t)[C_N] ...       v_n(t)[C_n] ...      v_1(t)[C_1]
%                   |                     |                    |
%   -            -  |                  -  |                 -  |
%   o---------------o---...---------------o---...--------------'
%
%
%   Inputs:
%       R: vector of resistances.
%       C: vector of capacitances.
%       ascending: name-value pair (default: false) to indicate if (R, C) 
%                  are in ascending (true; measurement across net 1) or 
%                  descending (false; measurement across net N) order.
%
%   Outputs:
%       If 1 output:
%           sys: state-space model for the RC ladder system, where sys.A is 
%                in tridiagonal form, sys.B and sys.C have one nonzero 
%                element.
%       If 4 outputs:
%           A: state-space A matrix.
%           B: state-space B matrix.
%           C: state-space C matrix.
%           D: state-space D matrix.
%           
%   See also:
%
%   $Author: BH$    $Date: 2022-04-10$  $Revision: 2$
%
%   REVISION 1: 2023-05-17
%       Added name-value parameter to produce models in ascending or
%       descending order.
%
%   REVISION 2: 2023-09-28
%       Allow output of LTI object or individual A, B, C, D matrices (to
%       enable compatibility with SDPVAR objects).
%
%   ©2023 ETH Zurich, Brett Hannigan; D-HEST; Biomedical and Mobile Health Technology (BMHT) Lab; Carlo Menon

    p = inputParser;
    addParameter(p, 'ascending', false, @islogical);
    parse(p, varargin{:});

    N = length(C);
    A = zeros(N);

    if ~iscolumn(R)
        R = R';
    end

    if p.Results.ascending
        A = diag(1./C)*(diag(-1./R) + diag([-1./R(2:end); 0]) + diag(1./R(2:end), 1) + diag(1./R(2:end), -1));
        B = [1/(R(1)*C(1)); zeros(N-1, 1)];
        C = [-1/R(1) zeros(1, N-1)];
        D = 1/R(1);
    else
    %    % Old method
    %     Crow_inv = diag(1./C);
    %     Rtri = tril(repmat(R, 1, N));
    %     N_inv = -eye(N) + diag(ones(N-1, 1), 1);
    % 
    %     A = Crow_inv/Rtri*N_inv;
    %     B = [zeros(N-1, 1); 1/(R(end)*C(end))];
    %     C = [zeros(1, N-1) -1/R(end)];
    %     D = 1/R(end);
        A = diag(1./C)*(diag(-1./R) + diag([0; -1./R(1:end-1)]) + diag(1./R(1:end-1), 1) + diag(1./R(1:end-1), -1));
        B = [zeros(N-1, 1); 1/(R(end)*C(end))];
        C = [zeros(1, N-1) -1/R(end)];
        D = 1/R(end);
    end

    if nargout==1
        varargout{1} = ss(A, B, C, D);
    else
        varargout{1} = A;
        varargout{2} = B;
        varargout{3} = C;
        varargout{4} = D;
    end

end