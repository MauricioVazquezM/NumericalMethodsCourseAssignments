%{

     Single Gaussian elimination step function

%}

%{
    
    COMMENTS: 

        - Gaussian elimination, also known as row reduction, is an 
        algorithm for solving systems of linear equations. It consists of 
        a sequence of operations performed on the corresponding matrix of 
        coefficients.

%}

%{

    INPUTS & OUTPUTS:

        - Input: An augmented matrix A ∈ R^(k,k+1) with A11 ̸= 0.

        - Output: A matrix B that has the same size as A and corresponds to a 
          single Gaussian elimination step applied to A using A11 as the 
          pivot (i.e. Bi1 = 0 for 2 ≤ i ≤ k).

%}


function B = elim(A)

    % Obtaining k-number of rows and m-number of columns
    [k, m] = size(A);

    % Checking that A11 is not zero
    if A(1,1) == 0
        
        % Pointing out the error
        error('The pivot element -A11- cannot be zero.');

    end

    % Checking that it's an augmented matrix
    if m ~= k + 1

        % Pointing out the error
        error('The input matrix is not augmented.');

    end

    % Defining B as A
    B = A; 

    % Starting looping in the second row, because the alrogithm is trying
    % to make all the elements below the pivot in the first column equal to
    % zero
    for i = 2:k

        % Defining the "factor" that will be used to eliminate the element 
        % in the first column of i row.
        fact = A(i,1) / A(1,1);

        % Row operation: the algorithm is subtracting 'factor' times the 
        % entire first row from the current 'ith' row to produce the 
        % new 'ith' row for matrix B
        B(i,:) = A(i,:) - fact * A(1,:);

    end

end
