%{

     My spline implementation algorithm

%}

%{

    INPUTS & OUTPUTS:

        - Input: x, y ∈ R^(n+1) and the elements of x are distinct.

        - Output: a, b, c, d ∈ R^n define the coefficients of n cubic 
        splines.

%}

function [a, b, c, d] = my_spline(x, y, dy)

    % Calculating the number of spline intervals, which is one less 
    % than the number of data points.
    n = length(x) - 1;

    % Computing the differences between consecutive x-values, which are the
    % widths of the spline 
    h = diff(x);
    
    % Initializing a square matrix A of size 4n by 4n and a vector B of 
    % size 4n. These will be used to set up a system of linear equations 
    % to solve for the spline coefficients
    A = zeros(4*n, 4*n);
    B = zeros(4*n, 1);
    
    % Setting up equations for the spline coefficients
    for i = 1:n

        % Setting up the coefficients for the spline function at x_i which 
        % should equal y_i
        A(i, 4*i-3:4*i) = [1, h(i), h(i)^2, h(i)^3]; % ( s(x_i) = y_i )

        % Setting up the coefficients for the spline function at x_{i+1} 
        % which should equal y_{i+1}
        A(n+i, 4*i-3:4*i) = [1, 0, 0, 0]; % ( s(x_{i+1}) = y_{i+1} )
        
        % Checking if the current interval is not the last one to set up
        % continuity conditions for inner intervals.
        if i < n

            % Setting up the conditions to ensure that the first derivative 
            % of the spline is continuous at the interior points
            A(2*n+i, 4*i-3:4*i) = [0, 1, 2*h(i), 3*h(i)^2];
            A(2*n+i, 4*(i+1)-3:4*(i+1)) = [0, -1, 0, 0];
            
            % Setting up the conditions to ensure that the second 
            % derivative of the spline is continuous at the interior 
            % points
            A(3*n+i-1, 4*i-3:4*i) = [0, 0, 2, 6*h(i)];
            A(3*n+i-1, 4*(i+1)-3:4*(i+1)) = [0, 0, -2, 0];

        end

    end
    
    % Setting up the conditions for the first derivative at the first and 
    % last points to match the given derivatives dy_0 and dy_n
    A(4*n-1, 1:4) = [0, 1, 0, 0]; % ( s'(x_0) = dy_0 )
    A(4*n, 4*n-3:4*n) = [0, 1, 2*h(n), 3*h(n)^2]; % ( s'(x_n) = dy_n )
    
    % Filling the B vector with the y-values and the specified derivatives 
    % at the endpoints
    B(1:n) = y(1:n);
    B(n+1:2*n) = y(2:n+1);
    B(4*n-1) = dy(1);
    B(4*n) = dy(n+1);
    
    % Solving the system of linear equations Ax = B for x, which contains 
    % the coefficients of the cubic splines in a vector
    coeffs = A\B;
    
    % Extracting the individual coefficients for a, b, c and d from 
    % the solution vector coeffs
    a = coeffs(1:4:end);
    b = coeffs(2:4:end);
    c = coeffs(3:4:end);
    d = coeffs(4:4:end);

end
