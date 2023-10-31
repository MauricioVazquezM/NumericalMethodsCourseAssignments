%{

     Gauss seidel method implementation

%}

%{

    INPUTS & OUTPUTS:

        - Input: The inputs are a matrix A ∈ R^(n,n) and a 
        vector b ∈ R^n.

        - Output: The output vector x ∈ R^ n is an approximate solution 
        of the linear equation Ax = b computed via Jacobi iteration niter 
        times. convergent is a Boolean value describing whether or not the
        iteration converged.

%}

function [x, niter, convergent] = gauss_seidel(A, b, tol)

    % Obtaining n-number of rows and m-number of column
    [n, m] = size(A);
    
    % Checking if A is square
    if m ~= n

        % Displaying error
        error('Matrix A is not square.');

    end
    
    % Checking dimensions of b
    if length(b) ~= n

        % Displaying error
        error('Dimension mismatch between A and b.');

    end

    % Initializing vector x solution
    x = zeros(n, 1);

    % Initializing number of iterations variable
    niter = 0;

    % Inializing convergence variable
    convergent = false;

    % Getting lower triangular matrix of A matrix
    L = tril(A);

    % Getting upper triangular matrix of A matrix
    U = triu(A, 1);

    while niter < 1000

        % Solving with our triupper solver implementation
        x_new = triu_solve(L, b - U*x);
        
        % Checking convergence
        if max(abs((x_new - x)./x_new)) < tol

            % Updating convergence variable
            convergent = true;

            break;

        end
        
        % Updating x solution vector
        x = x_new;

        % Updating niter variable
        niter = niter + 1;

    end

    % Checking if number of iterations has been reached
    if niter == 1000

        % Updating convergence variable
        convergent = false;

    end
    
end
