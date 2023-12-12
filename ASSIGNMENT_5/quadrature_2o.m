%{

     Quadrature function implementation

%}

function I = quadrature_2o(x, y)

%{

    INPUTS & OUTPUTS:

        - Inputs:
            1) x is a vector containing the points: 

                           ( a, (a+b)/2, b )

            2) y is a vector containing the function values: 

                           ( f(a), f((a+b)/2), f(b) )

        - Output:
            1) Integras I is obtained by the quadrature formula from the 
            previous exercise

%}

%{

    BUILT-IN MATLAB FUNCTIONS USED:

        - NA

%}

    % Calculating step size h
    h = (x(3) - x(1)) / 2; 

    % Calculating the integral with the previous excercise formula obtained
    I = (h/3) * (y(1) + 4*y(2) + y(3)); 

end
