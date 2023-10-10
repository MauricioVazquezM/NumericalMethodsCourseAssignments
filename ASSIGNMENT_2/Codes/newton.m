%{
    Newton Method
%}

%{
COMMENTS:
    Newton's method sucess depends on the initial choice of x0 and may not 
    converge if the initial guess is far from the root. Must be careful 
    when calculating the derivative of the function, since errors in its 
    calculation can affect the convergence of the method. 
%}

function [x, iterates, residuals] = newton(f, df, x0, maxiter, tol)

    % Define my auxiliar x
    x = x0;

    % Initialize iterates matrix
    iterates = zeros(length(x0) , maxiter); 

    % Initialize residuals vector
    residuals = zeros(maxiter , 1); 

    % For loop for Newton iteration
    for niter = 1: maxiter

        % Updating auxiliar
        fx = feval(f, x);
        dfx = feval(df, x);
        diff = -(dfx\fx');

        % Newton iterarion
        x = x + diff;

        % Store the current iterate and the residual
        % k-th column of the iterates matrix 
        iterates(: , niter) = x'; 

        % k-th position of the residuals vector
        residuals(niter) = norm(feval(f, x)); 

        % Check convergence depending on the tolerance
        if residuals(niter) < tol

            % Displaying result
            % fprintf('Convergence after %d iterations.\n', niter);
            
            % Store values
            iterates = iterates(: , 1:niter);
            residuals = residuals(1:niter);

            % Finish, function has found the root
            return; 

        end

    end

    % Max iterations reached, displaying result
    % fprintf('No convergence after %d iterations.\n', maxiter);

end



