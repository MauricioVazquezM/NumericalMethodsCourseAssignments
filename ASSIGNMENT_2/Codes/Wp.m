%{
    Lambert function principal branch script
%}

%{
COMMENTS:
    The Lambert W function is a valuable mathematical tool for solving 
    equations that involve exponential and logarithmic terms.
%}

function w = Wp(x0)
    
    if x0 >= 0

        % By hypothesis (Wikipedia source provided), defining the 
        % function and its derivative
        f = @(x) x * exp(x) - x0;
        df = @(x) exp(x) + x*exp(x);
    
        % Maximunm iterations for Newton's method, based on previous
        % excercises
        maxiter = 20;
    
        % Basing on the last excercise, setting tol = eps
        tol = eps;
    
        % Using newton function to find the solution
        [w, ~, ~] = newton(f, df, x0, maxiter, tol);

        % Display estimated root
        % fprintf('Estimated root: w= %d .\n', w);

    else
        
        % Max iterations reached, displaying result
        fprintf(['By hypothesis, no possible to obtain Wp or W0 because %d ' ...
            'is not >= 0.\n'], x0);

    end

end