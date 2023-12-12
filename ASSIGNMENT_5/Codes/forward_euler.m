%{

     Forward Euler function implementation

%}

function [t, x] = forward_euler(f, TSPAN, x0, h)

%{

    INPUTS & OUTPUTS:

        - Inputs:
            1) f is a function handle for a vector field f : R × R^n → R^n

            2) TSPAN is a time interval of the form [0, T]

            3) x0 ∈ R^n is the initial condition

            4) h > 0 is size of each Euler step

        - Output:
            1) t ∈ R^m is a uniform partition of [0, T] into subintervals 
            of width h

            2) x ∈ R^(m,n) is the approximate solution. The jth row of this 
            array is the forward Euler approximation of x(t_j)

%}

%{

    BUILT-IN MATLAB FUNCTIONS USED:

        - 'floor': round toward negative infinity
        - 'isrow': determine if input is row vector

%}

    % Determining number of steps needed
    numSteps = floor((TSPAN(2) - TSPAN(1)) / h);
    
    % Initializing time vector
    t = linspace(TSPAN(1), TSPAN(1) + numSteps*h, numSteps + 1);
    
    % Checking if x0 is a row vector and convert to column if necessary
    if isrow(x0)
        
        % Transposing vector
        x0 = x0';

    end
    
    % Initializing solution x-matrix 
    x = zeros(length(x0), numSteps + 1);

    % Setting the initial condition
    x(:, 1) = x0;

    % Iterating over each time step
    for i = 1:numSteps

        % Calculating the next value and ensure it's a column vector
        next_val = x(:, i) + h * f(t(i), x(:, i));

        % Checking if next_val is a row vector
        if isrow(next_val)

            % Transposing vector
            next_val = next_val';

        end

        % Storing the next value on the matrix solution
        x(:, i+1) = next_val;

    end

    % Transposing x if the initial condition was a row vector
    if isrow(x0)

        % Transposing solution matrix
        x = x';

    end

end

