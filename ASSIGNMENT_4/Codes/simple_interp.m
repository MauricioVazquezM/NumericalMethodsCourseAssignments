%{

     Simple interpolation implementation algorithm

%}

%{

    INPUTS & OUTPUTS:

        - Input: x, y ∈ R^(n+1) and the elements of x are distinct.

        - Output: a ∈ R^(n+1) contains the coefficients of the interpolant 
        for this data

%}

%{

    BUILT-IN MATLAB FUNCTIONS USED:

        - 'vander': function that returns the Vandermonde Matrix such 
        that its columns are powers of the vector v.

%}

function a = simple_interp(x, y)

    % Checking that x and y are vectors of the same size
    if length(x) ~= length(y)

        % Displaying error
        error('Vectors x and y must have the same length.');

    end
    
    % Checking that elements of x are distinct
    if length(unique(x)) ~= length(x)

        % Displaying error
        error('Elements of x must be distinct.');

    end

    % Generating the Vandermonde matrix for vector x
    V = vander(x);

    % Solving the linear system V * a = y for coefficients a
    a = V \ y;

end
