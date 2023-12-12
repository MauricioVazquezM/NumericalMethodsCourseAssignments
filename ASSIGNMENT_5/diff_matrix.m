%{

     Differentiation matrix function implementation

%}

function D = diff_matrix(a, b, n)

%{

    INPUTS & OUTPUTS:

        - Input: a, b & n. Where a < b define the interval for which the 
        function data appears. 'n' represents the number of points used to 
        discretize the interval [a, b].

        - Output: D ∈ R^(n,n) is the matrix satisfying Dy = y′ where 
                y = (f(x1), . . . , f(xn)) 
                            and 
                y′ ≈ (f′(x1), . . . , f′(xn))
        is the approximation obtained by centered differences 
        (at the interior nodes).

%}

%{

    BUILT-IN MATLAB FUNCTIONS USED:

        - NA

%}

    % Generating uniformly spaced nodes in the interval [a, b]
    % x = linspace(a, b, n);
    
    % Calculating the distance between adjacent nodes
    h = (b - a) / (n - 1);
    
    % Initializing the differentiation matrix D
    D = zeros(n, n);
    
    % Centered differences for interior nodes
    % Iterating over the interior nodes excluding the first and the last
    % nodes
    for i = 2:n-1

        % Applying the centered difference formula
        D(i, i-1) = -1 / (2*h);
        D(i, i+1) = 1 / (2*h);

    end
    
    % Appyling the forward difference for the first node
    D(1, 1) = -3 / (2*h);
    D(1, 2) = 2 / h;
    D(1, 3) = -1 / (2*h);
    
    % Appyling the backward difference for the last node
    D(n, n-2) = 1 / (2*h);
    D(n, n-1) = -2 / h;
    D(n, n) = 3 / (2*h);

end
