% sin_taylor function script

% Question 3: Auxiliar script with modified function

function ers=sin_taylor2(x,n)
    ers = zeros(size(n));
    n_v = 1:n;
    for i = 1: length(n_v)
        k = n_v(i);
        y = 0;
        for j = 0:k
            aux = (-1)^j*x^(2*j+1)/factorial(2*j+1);
            y = y + aux;
        end
    ers(i) = abs(sin(x) - y);
    end
end