% Trisection algorithm

% QUESTION 8 Write a MATLAB function called trisect.m which implements 
% the algorithm according to the prototype: 

function root = trisection(f, a, b, epsilon, number_iterations)
    syms x
    % Evaluation at the extremes of the interval
    fa = subs(f, x, a);
    fb = subs(f, x, b);
    % Validation
    if fa*fb > 0
        error('The initial interval [a, b] does not include a root.')
    else
        % Lets start...
        i = 0;
        while abs(a-b) < epsilon && i < number_iterations 
            if fa < fb % f fa is less than fb, 
                % it means that the sign of f at point a 
                % is different from the sign of f at point b
                aux1 = (2*a + b) / 3; % Potential trisection point
                f_aux1 = subs(f, x, aux1); % Evaluation on trisec point
                % Depending on the signs of fa, ev_aux1, 
                % and ev_aux2, the interval [a, b] is updated accordingly
                if fa*f_aux1 < 0
                    % Interval [a, b] is updated
                    b = aux1;
                else
                    aux2 = (a + 2*b) / 3; % Potential trisection point
                    f_aux2 = subs(f, x, aux2); % Evaluation on trisec point
                    if f_aux2*f_aux1 < 0
                        % Interval [a, b] is updated
                        a = aux1;
                        b = aux2;
                    else
                        % Interval [a, b] is updated
                        a = aux2;
                    end
                end
            else % If fb is less than fa, the logic is similar but with 
                 % reversed roles for a and b.
                aux1 = (a + 2*b) / 3; % Potential trisection point
                f_aux1 = subs(f, x, aux1); % Evaluation on trisec point
                % Depending on the signs of fa, ev_aux1, 
                % and ev_aux2, the interval [a, b] is updated accordingly
                if fb*f_aux1 < 0
                    % Interval [a, b] is updated
                    a = aux1;
                else
                    aux2 = (2*a + b) / 3; % Potential trisection point
                    f_aux2 = subs(f, x, aux2); % Evaluation on trisec point
                    if f_aux1*f_aux2 < 0
                        % Interval [a, b] is updated
                        a = aux2;
                        b = aux1;
                    else
                        % Interval [a, b] is updated
                        b = aux2;
                    end
                end
            end
            % Increments iteration counter +1
            i = i + 1;
        end
        % The code computes root as the midpoint of the final 
        % interval [a, b]
        root = (a + b) / 2;
    end
end