%% ASSIGNMENT 5

%{

Student data: 

    Name: Mauricio Vazquez Moran

    Student number: 2821507

%}

%% QUESTION 1

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in diff_matrix.m script 
            b) Proof of the correct performance of the function

%}

%{ 

    Part a)

        Check diff_matrix.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Cleaning workspace
clear all; close all;

% Defining the function and its derivative
f = @(x) 0.1*x.^2 + x.*sin(x);
f_p = @(x) 0.2*x + sin(x) + x.*cos(x);

% Defining the given interval
a = 0;
b = 10;

% Defining aux variable for range of n values to test
n_val = [10, 15, 20, 35, 50];

% Plotting

% Iterating over our n test values
for i = 1:length(n_val)

    % Setting the current n
    n = n_val(i);

    % Generating linear spaced points between a and b
    x = linspace(a, b, n);

    % Calling our function implementation
    D = diff_matrix(a, b, n);

    % Calculating the aproximative derivative of f at each point
    f_approx = D * f(x)';
    
    % Setting up a subplot for each iteration
    subplot(length(n_val), 1, i);

    % Plotting the actual derivative
    plot(x, f_p(x), 'b', x, f_approx, 'r--');

    % Adding legend 
    legend('Actual Derivative', 'Approximated Derivative');

    % Adding tittle to each subplot based on current n value
    title(['n = ', num2str(n)]);

end

% Displaying conclusion
fprintf(['By plotting the results, we can visually assess how well the ' ...
    'approximation converges to the actual derivative as n increase.  This ' ...
    'is because smaller intervals (smaller h) lead to better approximations ' ...
    'of the derivative. On the other hand, larger n also means more ' ...
    'computational work and memory usage. Therefore, a very large n' ...
    ' might not always be practical. \n'])

%% QUESTION 2

%{ 

COMMENTS: 

    - I will leave the quatradure formula procedure as pdf file. This is 
      because the procedure was easier to do on paper and pen than in 
      digital format. 

    - I will leave the proof procedure as a pdf file. This is because the 
      procedure was easier to do on paper and pen than in digital format. 
    
%}

% Cleaning
clear all, close all;

%{ 

    Part a)

        Check Exercise3_Quad_Formula_part.pdf file is on the same folder as 
        main.m

%}

%{ 

    Part b)

        Check Exercise3_Proof_part.pdf file is on the same folder as main.m

%}

%% QUESTION 3

%{ 

COMMENTS: 

    I am going to divide this question into 2 parts:
                a) Algorithm implementation in quadrature_2o.m script 
                b) Proof of the correct performance of the function
               
        
%}

%{ 

    Part a)

        Check quadrature_2o.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Cleaning
clear all, close all;

% Defining the function
f = @(x) 0.1*x.^2 + x.*sin(x);

% Defining the interval [a, b]
a = 0;
b = 10;

% Defining the number of subintervals
n = 100; 

% Creating a uniform partition of the interval [a, b]
x = linspace(a, b, n+1); 

% Evaluating the function at these points
y = f(x);

% Initializing the integral approx
I_tot = 0;

% Performing the quadrature for each consecutive triple
for i = 1:2:length(x)-2

    % Using quadrature function implementation
    I_tot = I_tot + quadrature_2o(x(i:i+2), y(i:i+2));

end

% Displaying the total integral approximation
disp('Integral value using quadrature function implementation: ');
disp(I_tot);

% Plotting the function and the quadrature points
figure;
fplot(f, [a b]);
hold on;

% Plotting the quadrature points
plot(x, y, 'square'); 

% Adding tittle
title('Function plot with quadrature points');

% Adding labels
xlabel('x');
ylabel('f(x)');

% Adding legend
legend('f(x)', 'Quadrature Points');
hold off;

% Verifying the convergence numerically
I_act = integral(f, a, b);

% Displaying results using MATLAB's integral function
disp('Integral value using MATLABs integral function: ');
disp(I_act);

% Displaying conclussions
error = I_tot - I_act;
fprintf(['Error of integral value between Quatrature function implementation ' ...
    'and MATLABs integral function: %.10f\n'], error);
fprintf('\n');
disp('Conclussion: ');
disp(['As we can see, the error is very small. Which allows us to conclude ' ...
    'that quadrature function implemented is a very good numerical ' ...
    'approximation for an integral.']);

%% QUESTION 4

%{

    COMMENTS: 

        I am going to divide this question into 3 parts:
            a) Algorithm implementation in forward_euler.m script 
            b) Proof of the correct performance of the function
            c) Plotting and comparison with 'ode45'

%}

%{

    ADDITIONAL REFERENCES:

        - Braun, M. (1983). Differential equations and their applications. 
          In Applied mathematical sciences. 
          https://doi.org/10.1007/978-1-4684-9229-3

        - Wikipedia contributors. (2023, September 26). Mass-spring-damper 
        model. Wikipedia. 
        https://en.wikipedia.org/wiki/Mass-spring-damper_model

%}

%{ 

    Part a)

        Check diff_matrix.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

%{

    INTRODUCTORY COMMENT TO EXAMPLE USED:

        - The mass-spring-damper model consists of discrete mass nodes 
        distributed throughout an object and interconnected via a network 
        of springs and dampers. This model is well-suited for modelling 
        object with complex material properties such as nonlinearity and 
        viscoelasticity.

%}

%{

    MASS-SPRING SYSTEM:
        
        - Consider a mass-spring system with mass m and spring constant k.
        Let x1 be displacement and x2 be the velocity. Then, the equations
        are: 

                dx1/dt = x2

                dx2/dt = -(k/m)x1

        - For simplicity, I choose m = 1 and k = 1, then:

                dx1/dt = x2

                dx2/dt = -x1

%}

% Cleaning env
clear all, close all;

% Defining equation
f = @(t, x) [x(2); -x(1)];

% Definining initial condition
x0 = [1; 0];  % Starting at x1 = 1 (displacement), x2 = 0 (velocity)

% Defining Time span
TSPAN = [0 10];

%{ 

    Part c)

        Plotting and comparison with 'ode45'

%}

% Defining different h values for testing
h_values = [0.1, 0.05, 0.01]; 

% For-loop for testing purposes
for h = h_values

    % Solving with our Forward Euler implementation
    [t_eul, x_eul] = forward_euler(f, TSPAN, x0, h);

    % Transposing x-matrix solution for plotting purporses
    x_eul = x_eul';

    % Using ode45 for comparison
    [t_ode45, x_ode45] = ode45(f, TSPAN, x0);

    % Plotting
    figure;
    plot(t_eul, x_eul(:, 1), 'b-', t_ode45, x_ode45(:, 1), 'r--');

    % Adding legend
    legend('Forward Euler implemenation', 'ode45 MATLABs BUILT-IN function');

    % Adding tittle
    title(['Solution Comparison with h = ', num2str(h)]);

    % Adding labels
    xlabel('Time');
    ylabel('Displacement');

end

% Displaying conclusion
fprintf(['As h gets smaller, the Forward Euler solution should get closer ' ...
    'to the ode45 solution. However, due to the first-order accuracy of ' ...
    'the Forward Euler method, very small step sizes might be needed to ' ...
    'achieve close agreement, and the method may not be stable for all ' ...
    'types of ODEs.\n'])

%% QUESTION 5

%{

    COMMENTS: 

        I am going to divide this question into 2 parts: 
            a) Simulation part and plotting part
            b) Discussing results
            c) Conclussion

%}

%{

    ADDITIONAL REFERENCES:

        - Weisstein, Eric W. "Lotka-Volterra Equations." From MathWorld--A 
          Wolfram Web Resource. 
          https://mathworld.wolfram.com/Lotka-VolterraEquations.html

        - Wikipedia contributors. (2023, 8 november). Lotka–Volterra 
          equations. Wikipedia. 
          https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations


%}

%{ 

    Part a)

        Simulation part and plotting part

%}

%{

    GENERAL INFORMATION:

        - The Lotka–Volterra equations, also known as the Lotka–Volterra 
        predator–prey model, are a pair of first-order nonlinear 
        differential equations, frequently used to describe the dynamics of 
        biological systems in which two species interact, one as a predator 
        and the other as prey. 

%}

%{

    -IDENTIFYING PREDATOR AND PREY ON THE GIVEN SYSTEM:

        - 1) Prey (A): The term A in A_dot = A − AB suggests that 
        the population of A grows linearly in the absence of B. However, 
        the negative interaction term −AB indicates that the presence of 
        B has a limiting effect on the growth of A, which is characteristic 
        of prey.

        - 2) Predator (B): The equation B_dot = AB − alpha*B includes a 
        positive interaction term AB, implying that the population of B 
        benefits from the presence of A, characteristic of a predator. The 
        term −alpha*B represents natural decay or death of the predator 
        species in the absence of the prey.

%}

%{

    X-MATRIX Solution interpretation:

        - The first column of x (x(:, 1)) represents the population of prey 
          at each time point

        - The second column of x (x(:, 2)) represents the population of 
          predators at each time point

%}

% Cleaning env
clear all, close all;

% Definining alpha value
alpha = 0.5;  

% Definining the ODE system
f = @(t, x) [x(1) - x(1)*x(2); x(1)*x(2) - alpha*x(2)];

% Defining initial conditions
init_cond = [1 1; 0.5 1.5; 1 0.5; 0.1 0.1];

% Defining time span
TSPAN = [0 30];

% Defininf step size
h = 0.01;

% Iterating over each initial condition value
for i = 1:size(init_cond, 1)

    % Getting initial condition value
    x0 = init_cond(i, :);

    % Using our Forward Euler function implementation
    [t, x] = forward_euler(f, TSPAN, x0, h);

    % Transposing x-matrix solution for plotting purporses
    x = x';

    % Plotting
    figure;
    plot(t, x(:, 1), 'b-', t, x(:, 2), 'r--');

    % Adding legend
    legend('Prey (A)', 'Predator (B)');

    % Adding tittle
     title(['Predator-Prey Dynamics (alpha = ', num2str(alpha), ', IC: A = ', num2str(x0(1)), ', B = ', num2str(x0(2)), ')']);

    % Adding labels
    xlabel('Time');
    ylabel('Population');

    % Adding grid
    grid on;

end

%{ 

    Part b)

        Discussing results

%}

% Questions to discuss
% Question & Answer 1
fprintf('\n')
fprintf(['Can the predator hunt the prey to extinction? If so what happens ' ...
    'to the prey?  \n'])
fprintf(['Answer: If the predator population B becomes too large, it may ' ...
    'overhunt the prey A, leading to a decline in both populations due to ' ...
    'lack of food.\n'])
fprintf('\n')

% Question & Answer 2
fprintf('\n')
fprintf([' Can the predator starve to extinction? If so what happens to the ' ...
    'prey? \n'])
fprintf(['Answer: The predator can starve to extinction if the prey ' ...
    'population is too low to sustain them. Following the extinction or ' ...
    'significant decline of predators, the prey population is likely to ' ...
    'grow, potentially leading to overpopulation, since their primary ' ...
    'limiting factor (predation) is removed\n'])
fprintf('\n')

% Question & Answer 3 
fprintf('\n')
fprintf('Can the species ever coexist? \n')
fprintf(['Answer: Coexistence is possible if the predator-prey interaction ' ...
    'reaches a dynamic equilibrium where neither species outcompetes the ' ...
    'other significantly.\n'])
fprintf('\n')

%{ 

    Part c)

        Conclussion

%}

fprintf('\n')
fprintf('CONCLUSION: \n')
fprintf(['This model encapsulates the complex interplay and delicate ' ...
    'balance between predator and prey populations. It highlights the ' ...
    'importance of key parameters and initial conditions in determining ' ...
    'the stability and sustainability of ecological systems. Real-world ' ...
    'ecosystems arecmore complex, but this model provides fundamental ' ...
    'insights into predator-prey interactions and their potential ' ...
    'implications.\n'])
fprintf('\n')





