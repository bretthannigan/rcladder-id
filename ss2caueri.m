function varargout = ss2caueri(varargin)
%SS2CAUERI Convert State-Space Model to Cauer I Canonical Form
%   Converts a generic continuous-time state-space model to the Cauer I
%   canonical form (realized by integrators) using a state transformation
%   computed from the Routh Array following the procedure from [1].
%
%   Inputs:
%       If 1 input provided:
%           sys: LTI model of system to be transformed.
%       If multiple inputs provided:
%           A: state-space A state transition matrix of system to be 
%              transformed.
%           B: state-space B input-state matrix of system to be 
%              transformed.
%           C: state-space C state-output matrix of system to be 
%              transformed.
%           D: state-space D feedthrough matrix of system to be transformed 
%              (equal to D matrix of original system).
%   Outputs:
%       If 1-2 outputs requested:
%           sys_t: LTI state-space  model of transformed system.
%           T: similarity transformation (defined as P2 in [1]), where:
%              A_t = T*A*T^-1
%              B_t = T*B
%              C_t = C*T^-1
%              D_t = D
%       If 3-4 outputs requested:
%           A_t: transformed state-space A state transition tridiagonal 
%                matrix.
%           B_t: transformed state-space B input-state matrix. 
%           C_t: transformed state-space C state-output matrix.
%           T: similarity transformation as described above.
%       If 5 outputs requested:
%           D_t: transformed state-space feedthrough matrix (equal to D
%                matrix of original system) is provided as the 4th output.
%           T: similarity transformation matrix is provided as the 5th
%              output.
%           
%   See also:
%       [1] Hwang, C., Nakano, B., & Asahi, T. (1984). Transformation of 
%           state-space model to Cauer I and II CFE canonical forms. 
%           International Journal of Systems Science, 15(7), 797–804. 
%           https://doi.org/10.1080/00207728408926600
%
%   $Author: BH$    $Date: 2022-03-30$  $Revision: 0$

    %% Manage inputs

    if nargin==1
        sys = ss(varargin{1});
    elseif nargin>=3
        A = varargin{1};
        B = varargin{2};
        C = varargin{3};
        if nargin==4
            D = varargin{4};
        else
            D = zeros(size(B, 2), size(C, 1));
        end
        sys = ss(A, B, C, D);
    end
    n = size(sys.A, 1);

    %% Generation of Routh Array
    % This section generates the triangular Routh array $R$ described in 
    % Equation (4) from [1] as well as the partial quotients from the continued
    % fraction expansion, $h$, shown in Equation (9).

    R = zeros(2*n+1);
    R(1,:) = [1 zeros(1, 2*n)];
    R(2,1:2*n) = arrayfun(@(x) sys.C*mpower(sys.A, x)*sys.B, [0:(2*n-1)]); % Markov parameters of the system.
    for i_row=3:(2*n+1)
        for i_col=1:(2*n+2-i_row)
            R(i_row,i_col) = R(i_row-2,i_col+1)-R(i_row-2,1)/R(i_row-1,1)*R(i_row-1,i_col+1);
        end
    end
    h = R(1:2*n,1)./R(2:(2*n+1),1);

    %% Transformation Matrix for Cauer Form 1 (Using Integrators)

    P = zeros(n, n); % P_20^T = 0
    P(1,:) = h(1)*sys.C; % Eq. (37)
    for i_P=2:size(P, 1) % The following lines are for Eq. (38)
        P(i_P,:) = h(2*i_P-2)*h(2*i_P-1)*P(i_P-1,:)*sys.A;
        if i_P<=2
            P(i_P,:) = P(i_P,:) + h(2*i_P-1)/h(2*i_P-3)*P(i_P-1,:);
        else
            P(i_P,:) = P(i_P,:) + ...
                       h(2*i_P-1)/h(2*i_P-3)*(1 + h(2*i_P-2)/h(2*i_P-4))*P(i_P-1,:) - ...
                       h(2*i_P-1)/h(2*i_P-5)*h(2*i_P-2)/h(2*i_P-4)*P(i_P-2,:);
        end
    end

    %% Manage outputs

    if nargout<3
        varargout{1} = ss(P*sys.A/P, P*sys.B, sys.C/P, sys.D, sys.Ts);
        if nargout==2
            varargout{2} = P;
        end
    elseif nargout>=3
        varargout{1} = P*sys.A/P;
        varargout{2} = P*sys.B;
        varargout{3} = sys.C/P;
        if nargout==4
            varargout{4} = P;
        elseif nargout==5
            varargout{4} = sys.D;
            varargout{5} = P;
        end
    end

end