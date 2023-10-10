%% ASSIGNMENT 2

%{

Student data: 
    Name: Mauricio Vazquez Moran
    Student number: 2821507
%}

%% QUESTION 1

%{ 
COMMENTS: 
    Check newton.m script 
%}

%% QUESTION 2

%{ 
COMMENTS: 
    - I used the common MATLABS plot for this excercise. This allow me to 
    set of coordinates connected by line segments, specify X and Y 
    as vectors of the same length
%}

% Cleaning
clear all, close all;

% Define the function to aproximate a solution
f = @(x)cosh(x)+cos(x)-3;

% Define the derivative of the function
df = @(x)sinh(x)-sin(x);

% Define x0 or initial guess
x0 = 3;

% Define number of maximum iterations
maxiter = 20;

% Define error tolerance
tol = eps;

% Calling newton method function
[x, iterates, residuals] = newton(f, df, x0, maxiter, tol);

% Calculating order of convergence
order_of_convergence = log(residuals(3:end) - x ./ residuals(2:end-1) - x ) ./ log(residuals(2:end-1) - x ./ residuals(1:end-2) - x);

%{
ORDER OF CONVERGENCE:
     As result we get vector containing the estimated order of convergence 
     for each iteration after the first two iterations.
%}

% Plot residuals
figure;
plot(1:length(residuals), residuals, 'square-', 'LineWidth', 1);

% Adding tittle
title('Residuals vs. Number of iterations');

% Adding labels
xlabel('Number of iterations');
ylabel('Residuals');

% Adding grid
grid on;

% Displaying results
fprintf('Estimated root: x = %.6f\n', x);
fprintf('Order of Convergence: %.4f\n', order_of_convergence(end));

% Displaying additional comments
fprintf(['Additional comments: It can be seen in the graph that the magnitude ' ...
    'of the residuals, in this exercise, falls rapidly as the iterations ' ...
    'increase. Thus, by the fifth iteration, the residual has practically ' ...
    'converged to zero. \n'])

%% QUESTION 3

%{ 
COMMENTS: 
    - I used the common MATLABS plot for this excercise. This allow me to 
    set of coordinates connected by line segments, specify X and Y 
    as vectors of the same length
%}

% Cleaning
clear all, close all;

% Define the function to aproximate a solution
f = @(x) cosh(x) + cos(x) - 2;

% Define the derivative of the function
df = @(x) sinh(x) - sin(x);

% Define x0 or initial guess
x0 = 3;

% Define number of maximum iterations
maxiter = 20;

% Define error tolerance
tol = eps;

% Calling newton method function
[x, iterates, residuals] = newton(f, df, x0, maxiter, tol);

% Calculate order of convergence
order_of_convergence = log(residuals(3:end) - x ./ residuals(2:end-1) - x ) ./ log(residuals(2:end-1) - x ./ residuals(1:end-2) - x);

%{
ORDER OF CONVERGENCE:
     As result we get vector containing the estimated order of convergence 
     for each iteration after the first two iterations.
%}

% Plot residuals
figure;
plot(1:length(residuals), residuals, 'square-', 'LineWidth', 1);

% Adding tittle
title('Residuals vs. Number of iterations');

% Adding labels
xlabel('Number of iterations');
ylabel('Residuals');

% Adding grid
grid on;

% Displaying results
fprintf('Estimated root: x = %.6f\n', x);
fprintf('Order of Convergence: %.4f\n', order_of_convergence(end));

% Displaying additional comments
fprintf(['Additional comments: It can be seen in the graph that the magnitude ' ...
    'of the residuals, in this exercise, falls rapidly as the iterations ' ...
    'increase. Thus, by the eighth iteration, the residual has practically ' ...
    'converged to zero. \n'])

%% QUESTION 4

%{
COMMENTS: 
    Check Wp.m script and Wm.m script 
%}

%% QUESTION 5

%{
    COMMENTS: 
        I am going to divide this section into 3 parts:
            1) Show Wp error 
            2) Show Wm error
            3) Plot the functions
%}

%{
    Part 1: Show Wp error
        By hypothesis (Wikipedia link provided) we get:
            y = Wp(x) if x >= 0.
        So to show the error, we would need to pass an argument that is
        magnitude is less tan zero.

        I used a try and catch block to attempt to call the functions 
        with the invalid input. 
%}

% Test Wp with an invalid input
% Defining an input smaller than zero
x_inv = - 0.1;

try

    % Evaluating in function
    r_inv = Wp(x_inv);
    
    % Displaying if its valid
    % fprintf('Wp(%f) = %f\n', x_inv, r_inv);

catch
    
    % Displaying if is invalid
    fprintf('Error!!! Wp(%f) is invalid.\n', x_inv);

end

%{
    Part 2: Show Wm error
        By hypothesis (Wikipedia link provided) we get:
            y = Wm(x) if -(1/e) =< x <0.
        So to show the error, we would need to pass an argument that is
        magnitude is ou of this interval.

        I used a try and catch block to attempt to call the functions 
        with the invalid input. 
%}

% Test Wp with an invalid input
% Defining an input smaller than zero
x_inv = - (1/exp(1)) - 0.1;

try

    % Evaluating in function
    r_inv = Wm(x_inv);
    
    % Displaying if its valid
    % fprintf('Wm(%f) = %f\n', x_inv, r_inv);

catch
    
    % Displaying if is invalid
    fprintf('Error!!! Wm(%f) is invalid.\n', x_inv);

end

%{
    Part 3: Plotting functions
        Using blue for the graph of Wp and magenta for the graph of Wm.
%}

% Cleaning
clear all, close all;

% Definining a range of x-values
x = linspace(-(1/exp(1)), 3, 1000);

% Initializing arrays to store y-values
y1 = zeros(size(x));
y2 = zeros(size(x));

% Calculate the corresponding y-values using the Lambert W functions
for i = 1:length(x)

    % By source provided, this is will be valuating in Wm
    if -(1/exp(1)) <= x(i) && x(i) < 0

        % Other branch time to evaluate (Wm)
        y2(i) = Wm(x(i));

        % Placeholder for values outside the domain of Wm
        y1(i) = NaN; 

    else
        
        % Principal branch time to evaluate (Wp)
        y1(i) = Wp(x(i));

        % Other branch (Wm)
        y2(i) = NaN; 

    end

end

% Creating the plot
figure;
plot(x, y1, 'b', 'LineWidth', 2);
hold on;
plot(x, y2, 'm', 'LineWidth', 2);

% Adding label names
xlabel('x');
ylabel('y');

% Adding title name
title('Lambert W Branches');
legend('Wp', 'Wm');
grid on;
hold off;

% Displaying additional comments
fprintf(['Additional comments: In the graph, you can clearly see the ' ...
    'branching of the Lambert function. On the one hand, the magenta color ' ...
    'allows us to see how the function comes from -1/e until it hits zero ' ...
    'through the non-main branch. After that, the main branch resumes ' ...
    '(blue color) for all values greater than zero. \n'])

%% QUESTION 6

%{
    Part a) Using MATLAB function 'roots':
    The roots function solves polynomial equations of the form:
        p1(x^n) +...+ pn(x) + pn+1=0. 
    Polynomial equations contain a single variable with nonnegative 
    exponents.
%}

% Defining coefficients on the necessary format the function needs
coef = [0.01, 1000, 0.01];

% Using the MATLAB function
roots_r1 = roots(coef);

% Selection the largest root approximation
r1 = max(roots_r1);

% Displaying results
fprintf('Estimated root for r1: x = %.15f\n', r1);

%{
    Part b) Using Newton Method:
    Newton's method sucess depends on the initial choice of x0 and may not 
    converge if the initial guess is far from the root.
%}

% Defining the given polynomial function and its derivative
f = @(x) 0.01*x.^2 + 1000*x + 0.01;
df = @(x) 0.02*x + 1000;

% Initial guess, given the first result
x0 = 0.0001; 

% Given by hypothesis, defining tolerance
tol = eps; 

% Given by some last excercise, defining iterations
maxiter = 20; 

% Call the newton function to find the root
[r2, ~, ~] = newton(f, df, x0, maxiter, tol);

% Displaying results
fprintf('Estimated root for r2: x = %.15f\n', r2);

%{
    Part c) Using quadratic formula:
    The quadratic formula helps us solve any quadratic equation. First, 
    we bring the equation to the form ax² + bx + c=0, where a, b, and c are 
    coefficients. Then, we plug these coefficients in the formula: 
        (-b±√(b²-4ac))/(2a) 
%}

% Defining coefficients
a = 0.01;
b = 1000;
c = 0.01;

% By formula, defining discriminant
discriminant = b^2 - 4*a*c;

% Obtaining root
r3 = (-b + sqrt(discriminant)) / (2*a);

% Displaying results
fprintf('Estimated root for r3: x = %.15f\n', r3);

%{
    Part d) Commenting results.
%}

% Displaying comments
fprintf(['Additional comments: First of all, in order to observe the change ' ...
    'between using one method or the other, we have to display the answer ' ...
    'to approximately 15 decimal places. By observing and analyzing the ' ...
    'result, we can see that both the use of the Matlab function "roots" ' ...
    'and the use of Newtons method have an identical magnitude when ' ...
    'displayed to 15 decimal places. This shows the accuracy of Newtons ' ...
    'algorithm for finding roots. On the other hand, we can see that when ' ...
    'using the quadratic formula, the result changes by several decimal ' ...
    'places. In other words, it is not as accurate as the first two methods ' ...
    'used in this exercise. \n']);

%% QUESTION 7

%{
    COMMENTS:
        - The constraint may be the equation that describes the boundary 
        of a region or it may not be.
        - Need to be careful with the fact that in some cases minimums and 
        maximums won't exist even though the method will seem to imply 
        that they do.
%}

%{
    EXCERCISE DIVISION:
        Part a) Finding minimum and maximum
        Part b) Checking graphically that its correct
%}

%{
    PART A) Finding minimum and maximum
%}

% Cleaning
clear all, close all;

% Defining the Lagrangian function
L = @(x) x(1).^3 + x(2).^3 + 3*x(1)*x(1) + x(3)*((x(1) - 3).^2 + (x(2) - 3).^2 - 9);

% Defining Lagrange multipliers system
f = @(x) [3*x(1).^2 + 3*x(2) + 2*x(3)*(x(1) - 3), 3*x(2).^2 + 3*x(1) + 2*x(3)*(x(2) - 3), (x(1) - 3).^2 + (x(2) - 3).^2 - 9];

% Defining equations derivatives
df = @(x) [6*x(1) + 2*x(3), 3, 2*(x(1) - 3); 3, 6*x(2) + 2*x(3), 2*(x(2) - 3); 2*(x(1) - 3), 2*(x(2) - 3), 0];

% Defining initial guess
init_guess = [0; 0; 0];

% Defining maximum number of iterations
maxiter = 100;

% Defining tolerance 
tol = eps;

% Newton function for the critical points
[x_crit, ~, ~] = newton(f, df, init_guess, maxiter, tol);

% Displaying critical points
fprintf('Critical point: [ x = %.5f , y = %.5f ] \n', x_crit(1), x_crit(2));

% Evaluating the objective function on critical point
xtrem_val = x_crit(1).^3 + x_crit(2).^3 + 3 * x_crit(1) * x_crit(2);

% Displaying extreme value
fprintf('Extreme Value: %.5f \n', xtrem_val)

%{
    PART B) Checking graphically that is correct
%}

% Defining plot
figure;

% Defining a range of x-values 
x1 = -4:0.4:4;

% Defining a range of y-values equal to x-values
y1 = x1;

% Meshgrid fos the matrix needed fot plotting
[X1,Y1] = meshgrid(x1);

% Calculating the corresponding z-values of F(x,y)
Z1 = X1.^3 + Y1.^3 + 3*X1.*Y1;

% Creating a 3D-plot of F
meshz(Z1);

%{
    Adding constraint figure
%}

% Defining a range of x-values 
x2 = 0:0.50:7;

% Defining a range of y-values equal to x-values
y2 = x2;

% Calculating the corresponding z-values of F(x,y)
z2 = (x2 - 3).^2 - (y2 - 3).^2 - 9;

% Adding a 3D-plot of the constraint
hold on;
plot3(x2, y2, z2,'-o', 'Color', 'r');

% Labeling the axes 
xlabel('x');
ylabel('y');
zlabel('z');

% Adding colorbar
colorbar;

% Adding a title
title('Plot of F(x, y) with constraint');

% Adding grig to the figure
grid on;

%{
    Adding critical point
%}

% Plot the critical point
scatter3(x_crit(1), x_crit(2), xtrem_val, 'Marker', 'o', 'Color', 'g', 'LineWidth', 2);

% Adding legend 
legend('F(x, y)', 'constraint', 'Critical point');

% Displaying comments
fprintf(['The results obtained and the results observed in the graph seem ' ...
    'correct. It was a great learning experience to observe how ' ...
    'using a method to find roots, such as Newtons method, of a certain ' ...
    'function, can be used to find critical points. One disadvantage, ' ...
    'which I realized in this excercise, is that a near local maxima or ' ...
    'local minima, due to oscillation, may cause to a slow ' ...
    'convergence.\n']);

%% QUESTION 8

%{
    COMMENTS:
        - Fast convergence: It converges fast, if it converges. Which 
        means, in most cases we get root (answer) in less number of steps.
        - It requires only one guess.
        - Derivation is more intuitive, which means it is easier to 
        understand its behaviour, when it is likely to converge and when 
        it is likely to diverge.
        
%}

%{
    EXCERCISE DIVISION:
        Part a) Calculating
        Part b) Plotting
%}

%{
    PART A) Calculating
%}

% Cleaning
clear all, close all;

% Defining the function g(x)
g = @(x) 4*x - 4*x.^2;

% Generating 1000 uniformly spaced in the interval [-0.1, 1]
interval = linspace(-0.1, 1, 1000);

% Defining number of maxmum iterations
max_iterations = 100;

% Defining tolerance 
tol = 1e-6;

% Initializing array to store roots 
con_roots = zeros(1000, 1);

% Defining the roots values
r1 = 0.45; 
r2 = 0.65; 
r3 = 0.3; 
r4 = 0.9; 

% Defining r^ array
roots = [r1, r2, r3, r4]; 

% Creating arrays to store legend entries
ent = cell(length(roots), 1);

% Creating arrays to store legend handles
hand = zeros(length(roots), 1);

% Iterating points
for i = 1:1000
    x0 = interval(i);
    
    % Applying Newton's method at each point
    [x, iterates, residuals] = newton(@(x)g(x) - x0, @(x) 4 - 8 * x, x0, max_iterations, tol);
    
    % Checking convergence
    if abs(g(x) - x0) < tol

        % Determining root that converged
        [~, index] = min(abs(x - roots));
        con_roots(i) = index;

    end

end

%{
    PART B) Plotting
%}

% Creating colormap array
colors = ['g', 'c', 'm', 'y'];

% Ploting 
for i = 1:1000
    if con_roots(i) > 0

        % Ploting each interval point
        plot(interval(i), 0, '.', 'Color', colors(con_roots(i)));
        hold on;

    end

end

% Customizing plot 
% Customizing x-label
xlabel('Interval [ -0.1, 1]');

% Customizing y-label
ylabel('Number of Iterations');

% Adding title
title(['Convergence of the 1000 uniformly spaced points on the interval ' ...
    '[-0.1, 1] with Newtons Method']);

% Updating legend entries and handles
for i = 1:length(roots)
    
    % Displaying legend entries
    ent{i} = sprintf('Root %d', i);

    % Plotting legend handles
    hand(i) = plot(NaN, NaN, 'o', 'Color', colors(i)); % Create empty plot for legend

end

% Creating the legend
legend(hand, ent);

% Putting grid
grid on;

% Displaying results
fprintf('Counts of points converging to each root:\n');

% Displaying for each root
for i = 1:length(roots)
    
    % Displaying each x amount of points converged to root x
    fprintf('Root %d: %d points\n', i, sum(con_roots == i));

end

% Displaying comments
fprintf(['Newtons method always converges to the root nearest the initial ' ...
    'condition: From the plot, we can observe that this claim is generally ' ...
    'true. Points tend to converge to the root nearest their initial ' ...
    'condition.\n']);
fprintf(['Each of the four roots attracts the same number of ' ...
    'iterations.\n']);
