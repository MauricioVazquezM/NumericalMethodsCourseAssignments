%{

     Solver for the equation triu(M)x = b

%}

%{

    INPUTS & OUTPUTS:

        - Input: The inputs are an invertible upper triangular 
        matrix A ∈ R^(n,n) and a vector b ∈ R^n. Otherwise raise an 
        exception.

        - Output: The output vector x ∈ R^n is an approximate solution of 
        the linear equation Ax = b computed via back substitution. You may 
        use the built-in dot function if you like.

%}

function x = triu_solve(A, b)

    % Determining the number of elements of the vector b
    n = length(b);

    % Checking if A is square
    if size(A,1) ~= size(A,2)

        % Displaying error
        error('Matrix A is not square.');

    end

    % Checking if A is upper triangular
    if any(tril(A,-1))

        % Displaying error
        error('Matrix A is not upper triangular.');

    end

    % Checking if A is invertible (diagonal elements must not be zero)
    if any(diag(A) == 0)

        % Displaying error
        error('Matrix A is not invertible (there is a zero diagonal element).');

    end

    % Initializating solution vector
    x = zeros(n, 1);

    % Performing back substitution process
    % starting from the last row and moving upwards
    for i = n:-1:1

        % Calculating the dot product
        x(i) = (b(i) - dot(A(i,i+1:n), x(i+1:n))) / A(i,i);

    end

end
