%{

     Linear solver implementation

%}

%{

    INPUTS & OUTPUTS:

        - Input: The inputs are a matrix A ∈ R^(n,n) that is upper 
        Hessenberg and 1 ≤ i < j ≤ n.

        - Output: The output G matrix is equal to the identity matrix 
        except the 2 × 2 submatrix of G consisting of only the i, j rows
        and columns which is a rotation matrix.

%}

%{

    BUILT-IN MATLAB FUNCTIONS USED:

        - 'eye': The eye function is a built-in MATLAB function that 
        generates identity matrices.

        - 'atan2': The atan2 function returns the four-quadrant inverse 
        tangent (tan^-1) of Y and X, which must be real.

%}

function G = givens_rotation(A, i, j)

    % Initializing G as an identity matrix of the same size as A
    G = eye(size(A));
    
    % Calculating the values for the Givens rotation
    % If A(j,i) is already 0, G remains the identity matrix
    if A(j,i) ~= 0

        % Calculating r (the sqrt(a^2 + b^2))
        thetha = atan2(A(j,i), A(i,i)); 

        % Calculating sine of the rotation angle
        s = sin(thetha);

        % Calculating cosine of the rotation angle
        c = cos(thetha);
        
        % Assigning the values to the appropriate elements in G
        G(i,i) = c;
        G(j,j) = c;
        G(i,j) = s;
        G(j,i) = -s;

    end
    
end
