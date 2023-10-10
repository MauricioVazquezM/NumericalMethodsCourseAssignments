%{
    Lambert function other branch script
%}

%{
COMMENT:
    Each branch corresponds to a different solution to the equation. The 
    principal branch (Wp) is commonly used, but other branches may be 
    needed in specific cases.
%}

function w = Wm(x0)
    
    if -(1/exp(1)) <= x0 && x0 < 0

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
        fprintf(['By hypothesis, no possible to obtain Wm because %d ' ...
            'is not  >= (-1/e) and < 0.\n'], x0);

    end

end