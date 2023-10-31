%{

     First interchange for the partial pivoting strategy function

%}

%{

    INPUTS & OUTPUTS:

        - Input: The input is any augmented matrix A ∈ R^(k,k+1) with 
        Ai1 ̸= 0 for at least one i ∈ {1, 2, . . . , k}.

        - Output: An matrix B that is identical to A except that exactly 
        two rows have been interchanged so that |B11| ≥ |Bi1| for 
        2 ≤ i ≤ k.

%}

function B = pivot(A)

    % Obtaining k-number of rows and m-number of columns
    [k, m] = size(A);

    % Checking that A is an augmented matrix
    if m ~= k + 1
        
        % Pointing out the error
        error('The input matrix is not augmented.');

    end
    
    % Checking that there is at least one non-zero element in the first 
    % column
    if all(A(:,1) == 0)

        % Pointing out the error
        error('All elements in the first column of matrix A are zero.');

    end

    % Defining B as A
    B = A; 

    % Finding the row index with the maximum absolute value in the first 
    % column
    [~, maxRowIndex] = max(abs(A(:,1)));

    % Swapping the rows of B-matrix to ensure that the element at top-left
    % position B11 has the largest absolute value among the elements of the
    % first column
    if maxRowIndex ~= 1

        % Storing the entire first row of matrix B
        temp = B(1,:);

        % Moving the row with the maximum absolute value to the first row
        B(1,:) = B(maxRowIndex,:);

        % Placing the original first row to the maxRowIndex Position
        B(maxRowIndex,:) = temp;

    end
    
end
