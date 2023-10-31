%{

     Jacobi implementation

%}

%{

    The Jacobi method is a method of solving a matrix equation on a matrix 
    that has no zeros along its main diagonal (Bronshtein and Semendyayev 
    1997, p. 892). Each diagonal element is solved for, and an approximate 
    value plugged in. The process is then iterated until it converges. This 
    algorithm is a stripped-down version of the Jacobi transformation 
    method of matrix diagonalization.

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

function [x, niter, convergent] = jacobi(A, b, tol)

    % Obtaining n-number of rows and m-number of column
    [n, m] = size(A);

    % Checking if A is square
    if m ~= n

        % Displaying error
        error('Matrix A is not square.');

    end

    % Checking if A and b are the same dimension
    if length(b) ~= n

        % Displaying error
        error('Dimension mismatch between A and b.');

    end

    % Initializating vector x with zeros
    x = zeros(n, 1);

    % Initializating auxiliar counter of iterations
    niter = 0;

    % Initializating variable convergent
    convergent = false;
    
    % Extracting dagonal matrix D
    D = diag(A);

    % Getting D inverse
    D_inv = 1./D;
    
    % Calculating R that is the 'rest' of A once its diagonal is removed
    R = A - diag(D);
    
    % Setting a max iteration to prevent infinite loop
    while niter < 1000 

        % Calculating x_new based on Jacobi iteration formula
        x_new = D_inv .* (b - (R * x));
        
        % Checking convergence based on relative difference
        if max(abs((x_new - x)./x_new)) < tol

            % Setting convergent in True value
            convergent = true;

            break;

        end

        % Updating x 
        x = x_new;

        % Updating iterattion auxiliar variable
        niter = niter + 1;

    end

    % Checking divergence (after max iterations)
    if niter == 1000

        % Setting convergent varible in false
        convergent = false;

    end

end
