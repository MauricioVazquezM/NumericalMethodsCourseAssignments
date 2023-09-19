% sin_taylor function 

% Question 1: Write a MATLAB function called sin_taylor.m 
% that approximates the function f(x) = sin(x) 
% using a Taylor series centered around 0

function y=sin_taylor(x,n)
    y=0;
    for i = 0:n
        aux = (-1)^i*x^(2*i+1)/factorial(2*i+1);
        y = y + aux;
    end
end


