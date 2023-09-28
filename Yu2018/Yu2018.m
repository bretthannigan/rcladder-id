function [sys_est, theta_est] = Yu2018(sys, A_theta, B_theta, C_theta, n_theta, varargin)
%YU2018 Identification of Structured State-Space Models
%   Identification of strucutred state-space models using Hankel matrix
%   rank constraint formulated as a difference of convex programming
%   problem.
%
%   Inputs:
%       sys: Black-box discrete-time system model.
%       A_theta: Function that takes parameter vector theta and returns the
%                structured state-space A matrix (affine in theta).
%       B_theta: Function that takes parameter vector theta and returns the
%                structured state-space B matrix  (affine in theta).
%       C_theta: Function that takes parameter vector theta and returns the
%                structured state-space C matrix (affine in theta).
%       n_theta: Number of unknown parameters, length of vector theta.
%       epsilon: (Optional, default: 0.1) name-value argument for the 
%                epsilon parameter, controlling the initialization error in 
%                Eq. (25) of [1].
%       rel_tol: (Optional, default: 1e-6) name-value argument for the
%                relative tolerance stopping criterion for the
%                difference-of-convex programming in Algorithm 1 of [1].
%       lambda: (Optional, default: 1e-3) name-value argument for the
%               penalty parameter (Lagrange multiplier) from Eq. (20) of
%               [1].
%       max_iter: (Optional, default: 100) name-value argument for the
%                 maximum number of iterations stopping criterion.
%       h: (Optional, default: order of sys) number of impulse response
%          samples forming the columns of the Hankel matrix (equivalently
%          the number of columns in the extended controllability matrix).
%       v: (Optional, default: order of sys) number of impulse response
%          samples forming the rows of the Hankel matrix (equivalently the
%          number of rows in the extended observability matrix).
%       verbose: (Optional, default: 0) Write to console the
%                optimization progress.
%
%   Outputs:
%       sys_est: Estimated strucutred system.
%       theta_est: Vector of estimated parameters. 
%           
%   See also:
%       [1] Yu, C., Ljung, L., & Verhaegen, M. (2018). Identification of 
%           structured state-space models. Automatica, 90, 54–61. 
%           https://doi.org/10.1016/j.automatica.2017.12.023
%
%   $Author: BH$    $Date: 2023-09-27$  $Revision: 0$

    n = size(sys.A, 1);

    p = inputParser();
    addParameter(p, 'epsilon', 0.1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'rel_tol', 1e-6, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'lambda', 1e-3, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'max_iter', 100, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'h', n, @(x) isinteger(x) && isscalar(x));
    addParameter(p, 'v', n, @(x) isinteger(x) && isscalar(x));
    addParameter(p, 'verbose', 1, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    
    v = p.Results.v;
    h = p.Results.h;
    
    nu = size(sys.B, 2);
    ny = size(sys.C, 1);
    np = n_theta;
    
    %% Problem Formulation
    Y = impulse2hankel(sys, v, h);
    
    X = sdpvar(ny*v, nu*h, 'full');
    O = sdpvar(ny*v, n, 'full');
    C = sdpvar(n, nu*h, 'full');
    
    Obar = sdpvar(ny*v, n, 'full');
    Cbar = sdpvar(n, nu*h, 'full');
    theta = sdpvar(np, 1);
    Atheta = A_theta(theta);
    Abar = sdpvar(n, n, 'full');
    
    constraints = [O(1:ny,:)==C_theta(theta)];
    constraints = [constraints, Obar(1:(v-1)*ny,:)==O((ny+1):v*ny,:)];
    constraints = [constraints, C(:,1:nu)==B_theta(theta)];
    constraints = [constraints, Cbar(:,1:(h-1)*nu)==C(:,(nu+1):h*nu)];
    Z = [X O Obar; C eye(n) Atheta; Cbar Atheta Abar];
    ym_options = sdpsettings('verbose', p.Results.verbose==2, 'solver', 'mosek');
    
    %% Initialization
    init_constraint = norm(Y - X, 'fro')^2<=p.Results.epsilon;
    init_obj = norm(Z, 'nuclear');
    
    init = optimize([constraints init_constraint], init_obj, ym_options);
    if init.problem~=0
        error('Initial solution failed.')
    end
    
    rel_change = @(theta, theta_prev) norm(theta-theta_prev, 2)/norm(theta_prev, 2);
    theta_prev = zeros(size(theta));
    
    %% Iterative Optimization
    i_iter = 1;
    while (rel_change(value(theta), theta_prev)>p.Results.rel_tol) && (i_iter <= p.Results.max_iter)
        if p.Results.verbose
            fprintf(print_progress(i_iter, p.Results.max_iter, rel_change(value(theta), theta_prev), p.Results.rel_tol, value(theta), mod(i_iter-1, 10)==0, true, false));
        end
        theta_prev = value(theta);
        [U, V] = svd_n(value(Z), n);
        obj = norm(Y - X, 'fro')^2 + p.Results.lambda*(norm(Z, 'nuclear') - trace(U'*Z*V));
        sol = optimize(constraints, obj, ym_options);
        if sol.problem~=0
            error('Iterative solution failed.')
        end
        i_iter = i_iter + 1;
    end
    if p.Results.verbose
        fprintf(print_progress(i_iter, p.Results.max_iter, rel_change(value(theta), theta_prev), p.Results.rel_tol, value(theta), mod(i_iter-1, 10)==0, true, true));
    end
    theta_est = value(theta);
    sys_est = ss(A_theta(theta_est), B_theta(theta_est), C_theta(theta_est), [], sys.Ts);
    
    %% Helper Functions
    function H = impulse2hankel(sys, v, h)
        % To construct the (v*ny) x (h*nu) block Hankel matrix from Eq. (11).
        M = impulse(sys, sys.Ts*(v+h-1));
        M = M(2:end,:,:);
        M = cellfun(@squeeze, num2cell(M, [2 3]), 'UniformOutput', false);
        H = cell2mat(M(hankel(1:v,v:v+h-1)));        
    end
    function Ob = obsv_ext(A, C, v) %#ok<DEFNU>
        % Modification of obsv.m to generate the v-steps extended
        % observability matrix from Eq. (13). Not used.
        n = size(A, 1);
        ny = size(C, 1);
        Ob = zeros(ny*v, n);
        Ob(1:ny,:) = C;
        for k=1:v-1
            Ob(k*ny+1:(k+1)*ny,:) = Ob((k-1)*ny+1:k*ny,:)*A;
        end
    end
    function Co = ctrb_ext(A, B, h) %#ok<DEFNU>
        % Modification of ctrb.m to generate the h-steps extended
        % controllability matrix from Eq. (13). Not used.
        n = size(A, 1);
        nu = size(B, 2);
        Co = zeros(n, n*nu);
        Co(:,1:nu) = B;
        for k=1:h-1
            Co(:,k*nu+1:(k+1)*nu) = A*Co(:,(k-1)*nu+1:k*nu);
        end
    end
    function [U, V] = svd_n(Z, n)
        [U, ~, V] = svd(Z, 'econ');
        U = U(:,1:n);
        V = V(:,1:n);
    end
    function str = print_progress(iter, max_iter, rel_change, rel_tol, theta, is_header, is_data, is_footer)
        str_theta = num2str(theta', '% 10.3e ');
        str_data = sprintf('\n| %5i | %8i |  %9.3e | %9.3e | %s |', iter, max_iter, rel_change, rel_tol, str_theta);
        str_horzline = repmat('-', 1, length(str_data)-1);
        str = '';
        if is_header
            str_columns = '| iter  | max_iter | rel_change | tol       | theta ';
            str = [newline(), str_horzline, newline(), sprintf(str_columns), repmat(' ', 1, length(str_data)-length(str_columns)-2), sprintf('|\n'), str_horzline];
        end
        if is_data
            str = [str, str_data];
        end
        if is_footer
            str = [str, newline(), str_horzline, newline()];
        end
    end
end