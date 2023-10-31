%{

     Linear solver implementation

%}

%{

    INPUTS & OUTPUTS:

        - Input: The inputs are an invertible matrix A ∈ R^(n,n) and a 
        vector b ∈ R^n.

        - Output: The output vector x ∈ R^n is a solution of the linear 
        equation Ax = b computed using Gaussian elimination.

%}

function x = lin_solve(A, b)

    % Obtaining n-number of rows and m-number of columns
    [n, m] = size(A);

    % Checking if Matrix A is square
    if m ~= n

        % Displaying error
        error('Matrix A is not square.');

    end

    % Augmenting the matrix A with vector b
    AugmentedMatrix = [A b];

    % Checking for singularity
    if any(diag(AugmentedMatrix))==0
        
        % Displaying error
        error('Matrix A is not invertible (no zero diagonal elements).');

    end

    % Gaussian elimination with partial pivoting
    for i = 1:n-1

        % Partial pivoting for column i
        AugmentedMatrix(i:end, i:end) = pivot(AugmentedMatrix(i:end, i:end));

        % Eliminating below the pivot
        AugmentedMatrix(i:end, i:end) = elim(AugmentedMatrix(i:end, i:end));

    end

    % Extracting the upper triangular matrix and modified b vector
    TM = AugmentedMatrix(:, 1:end-1);
    b_mod = AugmentedMatrix(:, end);

    % Back substitution using triu_solve
    x = triu_solve(TM, b_mod);

end


