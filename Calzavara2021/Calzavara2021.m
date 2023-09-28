function [G, sys_rc, T] = Calzavara2021(varargin)
%CALZAVARA2021 Structured Identification of RC Models
%   Given an RC network modelled by: 
%       G*x'(t) = S*x(t) + B*u(t)
%       y(t) = C*x(t)
%   with (n x 1) state vector x, (m x 1) input vector u, (p x 1) output
%   vector y, symmetric adjacency matrix S, and diagonal admittance matrix
%   G = diag([C_1, C_2, ..., C_n]). 
%
%   Inputs:
%       sys: identified state-space model with general form.
%
%   Outputs:
%       sys_rc: structured state-space model.
%           
%   See also:
%       [1] Calzavara, G., Consolini, L., & Kavaja, J. (2021). Structured
%           identification for network reconstruction of RC-models. Systems
%           and Control Letters, 147, 104849. 
%           https://doi.org/10.1016/j.sysconle.2020.104849
%
%   $Author: BH$    $Date: 2022-04-13$  $Revision: 0$

    p = inputParser;
    addRequired(p, 'sys', @(x) isa(x, 'lti'));
    addRequired(p, 'C', @ismatrix);
    addOptional(p, 'B', [], @(x) ismatrix(x) || isempty(x));
    parse(p, varargin{:});
    
    sys = p.Results.sys;
    Ahat = sys.A;
    Bhat = sys.B;
    Chat = sys.C; 
    Dhat = sys.D;
    B = p.Results.B;
    C = p.Results.C;
    
    nx = size(Ahat, 1);
    nu = size(Bhat, 2);
    ny = size(Chat, 1);
    if isempty(B)
        B = zeros(nx, 0);
    end
    constraints = [];

    [V, Lambda] = eig(Ahat);
    if length(uniquetol(diag(Lambda), 10*eps))==nx 
        % Eigenvalues of Ahat are unique, so D is diagonal.
        D = diag(sdpvar(nx, 1));
    else
        D = sdpvar(nx, nx);
    end
    H = diag(sdpvar(nx, 1));

    if size(Bhat, 2)==0%all(Bhat==0) % System is autonomous, can use Remark 3 of [1]
        constraints = [D>=0, Lambda*D==D*Lambda, H>=0, Chat*V*D*V'*Chat'==C*H*C'];
        ym_options = sdpsettings('verbose', true, 'solver', 'mosek');
    else
        W = inv(V);
        invD = sdpvar(D);
        constraints = [D>=0, Lambda*D==D*Lambda, H>=0, Chat*V*D*V'*Chat'==C*H*C', Bhat'*W'*invD*W*Bhat==B'*H*B, Chat*Bhat==C*H*B, D*invD==eye(nx)]
        ym_options = sdpsettings('verbose', true, 'solver', 'bmibnb');
    end
    
    solution = optimize(constraints, [], ym_options)
    if solution.problem==0
        H = value(H);
        G = inv(H);
        D = value(D);
        M = V*D*V';
        P = sqrtm(M);
        W = P*Chat';
        Z = sqrt(G)\C';
        if rank(Z)~=size(Z, 1) % Proposition 5 of [1].
            Q = W/Z + null(W')*null(Z')'; % Take a particular solution.
        else
            Q = W/Z;
        end
        T = P*Q*sqrt(G);
        sys_rc = ss2ss(sys, T);
    else
        solution.info
        yalmiperror(solution.problem)
    end
    
end