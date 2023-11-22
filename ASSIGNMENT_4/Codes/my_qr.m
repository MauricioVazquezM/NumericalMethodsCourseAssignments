%{

     My QR Algorithm implementation

%}

%{

    INPUTS & OUTPUTS:

        - Input: The inputs are a matrix A ∈ R^(n,n) that is upper 
        Hessenberg and 1 ≤ i < j ≤ n.

        - Output: The QR factorization.

%}

%{

    BUILT-IN MATLAB FUNCTIONS USED:

        - 'eye': The eye function is a built-in MATLAB function that 
        generates identity matrices.

%}

function [Q, R] = my_qr(A)

    % Getting size of A
    [m, n] = size(A);
    
    % Initializing Q as an identity matrix and R as A
    Q = eye(n);
    R = A;
    
    % Applying Givens rotations to zero out the elements below the main 
    % diagonal of R
    for j = 1:n-1

        for i = j+1:n

            if R(i,j) ~= 0

                % Calculating the Givens rotation matrix G that zeroes 
                % out R(i,j)
                G = givens_rotation(R, j, i);
                
                % Applying G to R from the left (pre-multiplication)
                R = G * R;
                
                % Accumulating the Givens rotations into Q
                Q = Q * G';

            end

        end

    end
    
end
