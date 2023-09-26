% sin_taylor function script

% Question 5: Auxiliar script with modified function.

function tay_aprox = sin_taylor3(x,n)
    y=0;
    tay_aprox = zeros(1, length(x));
    for i = 0:n
        aux = (-1)^i*x.^(2*i+1)/factorial(2*i+1);
        y = y + aux;
        tay_aprox = tay_aprox + aux;
    end
end