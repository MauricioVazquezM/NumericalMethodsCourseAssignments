%% ASSIGNMENT 3

%{

Student data: 

    Name: Mauricio Vazquez Moran

    Student number: 2821507

%}

%% QUESTION 1

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in elim.m script 
            b) Proof of the correct performance of the function

%}

%{ 

    Part a)

        Check elim.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Cleaning
clear all, close all;

% Defining the M-matrix 
M = [-1 11 -1 3; 0 3 -1 8; 10 -1 2 0; 2 -1 10 -1];

% Defining the b-vector
b = [25; 15; 6; -11];

% Augmenting M with b
AugmentedM = [M, b];

% Applying single Gaussian elimination step function
try

    % Using our implementation in elim.m script
    B = elim(AugmentedM);

    % Displaying result
    disp('Result after single Gaussian elimination step function:');
    disp(B);

catch exe

    % Displaying error
    disp(['Error: ', exe.message]);

end

% Testing exception for A11 = 0
TestM = [0 1 2 3 4; 5 6 7 8 9; 10 11 12 13 14; 15 16 17 18 19];

try
    
    % Using our implementation in elim.m script
    disp('Testing with matrix having A11 = 0:');
    B = elim(TestM);

catch exe

    % Displaying error
    disp(['Error: ', exe.message]);

end

% Displaying conclusion
fprintf(['Given the couple of tests done, it can be concluded that the ' ...
    'single Gaussian elimination step function works without any errors.\n'])

%% QUESTION 2

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in pivot.m script 
            b) Proof of the correct performance of the function

%}

%{ 

    Part a)

        Check pivot.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Cleaning
clear all, close all;

% Defining the M-matrix 
M = [-1 11 -1 3; 0 3 -1 8; 10 -1 2 0; 2 -1 10 -1];

% Defining the b-vector
b = [25; 15; 6; -11];

% Augmenting M with b
AugmentedM = [M, b];
% p = AugmentedM(2:end, 2:end);
% disp(p)

% Applying pivot function

try
    
    % Using our implementation in pivot.m script
    B_piv = pivot(AugmentedM);

    % Displaying the result
    disp('Result after partial pivoting function:');
    disp(B_piv);

catch exe

    % Displaying error
    disp(['Pivot Error: ', exe.message]);

end

% Testing exception for all zeros in the first column

% Defininf th matrix with all zeros in the first column
TestM0s = [0 1 2 3 4; 0 5 6 7 8; 0 9 10 11 12; 0 13 14 15 16];

try

    % Using our implementation in pivot.m script
    disp('Testing with matrix having all zeros in the first column:');
    B = pivot(TestM0s);

catch exe

    % Displaying error
    disp(['Pivot Error: ', exe.message]);

end

% Displaying conclusion
fprintf(['Given the couple of tests done, it can be concluded that the ' ...
    'partial pivoting function works without any errors.\n'])

%% QUESTION 3

%{ 
COMMENTS: 
    
%}

% Cleaning
clear all, close all;

% Defining the M-matrix 
M = [-1 11 -1 3; 0 3 -1 8; 10 -1 2 0; 2 -1 10 -1];

% Defining the b-vector
b = [25; 15; 6; -11];

% Testing using triu_solve implementation

try

    % Solving for the upper triangular part of M (trium(M))
    x_sol = triu_solve(triu(M), b);
    
    % Displaying solution
    disp('Solution x:');
    disp(x_sol);
    
    % Computing the residual norm
    res = b - triu(M) * x_sol;

    % Calculating the norm of the residual vector
    norm_res = norm(res);

    % Displaying results
    disp('Residual norm:');
    disp(norm_res);

catch exe
    
    % Displaying error
    disp(['Triu_solvers Error: ', exe.message]);

end


% Testing exception for non-square matrix

% Defining a non square matrix
TmatNonSq = [1 2; 0 3; 0 0];

% Defining our test b-vector
TVector = [3; 4; 5];

try

    % Trying to solve
    disp('Testing with non-square matrix:');
    x_test = triu_solve(TmatNonSq, TVector);

catch exe

    % Displaying error
    disp(['Triu_solver Error: ', exe.message]);

end

% Displaying conclusion
fprintf(['It can be concluded that the solver ' ...
    'implemented for triangular upper matrix works correctly, given that the' ...
    'residual norm result is 0.\n'])

%% QUESTION 4

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in line_solve.m script 
            b) Proof of the correct performance of the function

%}

%{ 

    Part a)

        Check line_solve.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Cleaning
clear all, close all;

% Defining the matrix M
M = [-1, 11, -1, 3; 0, 3, -1, 8; 10, -1, 2, 0;2, -1, 10, -1];

% Defining the b vector
b = [25; 15; 6; -11];

% Solving the equation Mx=b
try
    
    % Solving using our linear solver implementation
    x_sol = lin_solve(M, b);

    % Displaying results
    disp('Solution using linear solver implementation:');
    disp(x_sol);

    % Computing the residual norm
    res = b - M * x_sol;
    norm_res = norm(res);

    % Displaying results of residual norm
    disp('Residual norm using linear solver implementation:');
    disp(norm_res);

catch exe

    % Displaying error
    disp(['Linear solver Error: ', exe.message]);

end

% Testing  for a singular matrix

% Defining the singular matrix
SingularMatrix = [1, 1; 1, 1];

% Defining the vector
TestVector = [2; 2];

try

    % Trying to solve
    disp('Testing with singular matrix:');
    x_test = lin_solve(SingularMatrix, TestVector);

catch exe

    % Displaying error
    disp(['Linear solver Error: ', exe.message]);

end

% Displaying conclusion
fprintf(['It can be concluded that the linear solver implemented, using as ' ...
    'auxiliary functions those already implemented in the previous ' ...
    'exercises, worked correctly given that the result of the residual ' ...
    'norm is zero..\n'])

%% QUESTION 5

%{

    COMMENTS: 

        I am going to divide this question into 4 parts:
            a) Algorithm implementation in lin_solve_nopivot.m script 
            b) Main script work part
            c) Plotting results
            d) Discussing results

%}

%{ 

    Part a)

        Check lin_solve_npivot.m script 

%}

%{ 

    Part b)

        Main script work part

%}

% Cleaning
clear all, close all;

% Initializing size parameter
n = 100;

% Initializing number of linear problems parameter
num_lp = 1000;

% Initializing vectors to store relative residual errors
rel_res_err_lin_solve = zeros(num_lp, 1);
rel_res_err_nopivot = zeros(num_lp, 1);
% rel_res_err_MATLAB = zeros(num_lp, 1);

% Initializing vectors to store relative errors
rel_err_lin_solve = zeros(num_lp, 1);
rel_err_nopivot = zeros(num_lp, 1);
% rel_err_MATLAB = zeros(num_lp, 1);


for i = 1:num_lp

    % Generating random matrix A of nxn size
    A = randn(n, n);

    % While-loop checking the rank of A. If A is singular(rank less than n)
    % a new random matrix is generated until get a non-singular matrix(full
    % rank matrix).
    while rank(A) < n

        A = randn(n, n);

    end

    % Generating a random vector b
    b = randn(n, 1);

    % Checking the matrix for answer questions
    %AugmentedM = [A, b];
    %disp(AugmentedM)

    % Using linear solver implementation
    x1 = lin_solve(A, b);

    % Using linear solver without partial pivoting implementation
    x2 = lin_solve_nopivot(A, b);

    % Using MATLAB's backslash operator
    x_ML = A\b;

    % Computing relative residual errors
    rel_res_err_lin_solve(i) = norm(A*x1 - b, inf) / norm(b, inf);
    rel_res_err_nopivot(i) = norm(A*x2 - b, inf) / norm(b, inf);
    % rel_res_err_MATLAB (i) = norm(A*x_ML - b,inf) / norm(b, inf);

    % Computing residual errors
    rel_err_lin_solve(i) = norm(x_ML - x1, inf) / norm(x_ML, inf);
    rel_err_nopivot(i) = norm(x_ML - x2, inf) / norm(x_ML, inf);
    % rel_err_MATLAB (i) = norm(x_ML - x_ML,inf) / norm(x_ML, inf);

end

%{ 

    Part c)

        Plotting results 

%}

% Plotting the relative residual errors
figure;

% Adding values to graph
plot(rel_res_err_lin_solve,'-o','DisplayName','Relative Residual error with Partial Pivoting');
hold on;
plot(rel_res_err_nopivot,'-x','DisplayName','Relative Residual error without Partial Pivoting');
% plot(rel_res_err_MATLAB,'-s','DisplayName','MATLAB backlash Relative Residual error');

% Adding title
title('Relative Residual Errors');
% Adding x label name
xlabel('Number of trials');
% Adding y label name
ylabel('Error value');

% Adding legend
legend;
hold off;

% Plotting the relative errors
figure;

% Adding values to graph
plot(rel_err_lin_solve,'-o','DisplayName','Relative error with Partial Pivoting');
hold on;
plot(rel_err_nopivot,'-x','DisplayName','Relative error without Partial Pivoting');
% plot(rel_err_MATLAB,'-s','DisplayName','MATLAB backlash Relative error');

% Adding title
title('Relative Errors');
% Adding x label name
xlabel('Number of trials');
% Adding y label name
ylabel('Error value');

% Adding legend
legend;
hold off;

%{ 

    Part d)

        Discussing results

%}

% Questions to discuss

% Question & Answer 1
fprintf('\n')
fprintf('Does Gaussian elimination with/without pivoting reliably produce small residuals? \n')
fprintf(['Answer: Gaussian elimination with pivoting is more reliable in producing ' ...
    'small residuals than without pivoting. The pivoting strategy helps to handle numerical ' ...
    'inaccuracies that arise when manipulating numbers, especially when there are very small or zero ' ...
    'diagonal elements. Without pivoting, if theres a zero or a very small number in the diagonal, this ' ...
    'can cause large errors and hence large residuals. \n'])
fprintf('\n')

% Question & Answer 2
fprintf('\n')
fprintf('Does Gaussian elimination with/without pivoting reliably produce accurate approximations? \n')
fprintf(['Answer: Gaussian elimination with pivoting is usually more reliable in producing accurate ' ...
    'approximations. Pivoting ensures that the algorithm doesnt get derailed ' ...
    'by numerical inaccuracies. Without pivoting, the algorithm is susceptible to errors, especially ' ...
    'in cases where the matrix has nearly linearly dependent rows. \n'])
fprintf('\n')

% Question & Answer 3 
fprintf('\n')
fprintf('How often does pivoting improve each of the errors? \n')
fprintf(['Answer:  Given the test done, pivoting tends to improve errors more frequently when the matrix elements vary ' ...
    'widely in magnitude or when the matrix is nearly singular. \n'])
fprintf('\n')

% Question & Answer 4
fprintf('\n')
fprintf('When it does improve one or both errors, how large is that improvement? \n')
fprintf(['Answer: Given the tests done, depends on the specific system of equations. For some systems, the ' ...
    'improvement might be minimal, but for others (especially ill-conditioned systems), the ' ...
    'improvement can be significant. \n'])
fprintf('\n')

%% QUESTION 6

%{
    
    COMMENTS:

        I will leave the proof as a comment in this block.
%}

% Cleaning
clear all, close all;

%{
    PROOF:

    Let's consider a vector 'b' in the real numbers and a square matrix 'A' 
    of size nxn. 
    Assume that the diagonal elements of 'A' are strictly greater in 
    magnitude than the sum of the magnitudes of the other elements in their 
    respective rows.
    
    Based on the Gershgorin circle theorem, this ensures that 'A' is 
    invertible. This means that the equation Ax = b has a unique solution, 
    which we'll call 'x_hat'.
    To demonstrate that the sequence produced by the Jacobi iteration 
    converges, we'll refer to a theorem discussed in class about linear 
    iterations. 
    This theorem states that if the infinity norm of the matrix raised to 
    some power 'p' approaches zero, then the matrix itself converges.
    
    For the Jacobi method, the matrix in question is 'M', defined as the 
    inverse of the diagonal of 'A' times the sum of the lower and upper 
    triangular matrices of 'A'. In other words, M = D^(-1) * (L + U), where 
    'D' is the diagonal matrix of 'A', and 'L' and 'U' are its lower and 
    upper triangular matrices, respectively.
    
    The infinity norm of 'M' is the maximum row sum of its absolute values. 
    Using our assumption about 'A', we can show that this norm is strictly 
    less than 1.
    
    Since the infinity norm of 'M' is less than 1, as we raise 'M' to 
    higher and higher powers, its norm approaches zero. This means that 'M' 
    is convergent, and so is the sequence produced by the Jacobi iteration.
    
    This concludes the proof.

%}

%% QUESTION 7

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in jacobi.m script 
            b) Proof of the correct performance of the function

%}

%{ 

    Part a)

        Check jacobi.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Testing for diagonally dominant matrix

% Defining a diagonally dominant matrix A
A = [4, -1, 0, 0; -1, 4, -1, 0; 0, -1, 4, -1; 0, 0, -1, 4];

% Defining our b vector
b = [15; 10; 10; 10];

% Using the Jacobi method implementation
[x, niter, convergent] = jacobi(A, b, 1e-6);

% Displaying results
disp('Solution for the diagonally dominant matrix:');
disp(x);

% Displaying number of iterations
disp(['Iterations: ', num2str(niter)]);

% Displaying convergence
disp(['Convergent (1 if True, 0 if false): ', num2str(convergent)]);


% Testing on a non diagonally dominant matrix 

% Defining M matrix
M = [-1, 11, -1, 3; 0, 3, -1, 8; 10, -1, 2, 0; 2, -1, 10, -1];

% Defining b-vector
b_M = [25; 15; 6; -11];

% Using the Jacobi method implementation
[x_M, niter_M, convergent_M] = jacobi(M, b_M, 1e-6);

% Displaying solution
disp('Solution for matrix M:');
disp(x_M);

% Displaying number of iterartions
disp(['Iterations: ', num2str(niter_M)]);

% Displaying convergence
disp(['Convergent (1 if True, 0 if false): ', num2str(convergent_M)]);

% Displaying conclusion
fprintf(['The expected outcome is that the Jacobi iteration will converge ' ...
    'quickly for the diagonally dominant matrix A but might not converge or ' ...
    'may take longer for the non-diagonally dominant matrix M. The displayed ' ...
    'results help to understand and compare the performance of the Jacobi ' ...
    'method on these two matrices.\n'])

%% QUESTION 8

%{

    COMMENTS: 

        I am going to divide this question into 4 parts:
            a) Algorithm implementation in gauss_seidel.m script 
            b) Main script work part
            c) Plotting results
            d) Answering questions
            e) Repeating the experiment

%}

%{ 

    Part a)

        Algorithm implementation in gauss_seidel.m script

%}

%{ 

    Part b)

        Main script work part

%}

% Cleaning
clear all, close all;

% Initializing n an number of trials variables
n = 100;
num_trials = 1000;

% Initializing that will store the number of iterations each method took 
% for each trial and whether or not they converged.
iter_jacobi = zeros(num_trials, 1);
iter_gauss_seidel = zeros(num_trials, 1);
conv_jacobi = zeros(num_trials, 1);
conv_gauss_seidel = zeros(num_trials, 1);

% For-loop for trials
for i = 1:num_trials

    % Generating a random matrix A
    A = randn(n, n);

    % While-loop checking the rank of A. If A is singular(rank less than n)
    % a new random matrix is generated until get a non-singular matrix(full
    % rank matrix).
    while rank(A) < n

        A = randn(n, n);

    end

    % Generating a random vector b
    b = randn(n, 1);

    % Calling the jacobi and gauss_seidel functions implementations to 
    % solve the linear system
    [~, niter_jacobi, convergent_jacobi] = jacobi(A, b, 1e-6);
    [~, niter_gauss_seidel, convergent_gauss_seidel] = gauss_seidel(A, b, 1e-6);

    % Storing the results of each trial in the arrays initialized 
    iter_jacobi(i) = niter_jacobi;
    iter_gauss_seidel(i) = niter_gauss_seidel;
    conv_jacobi(i) = convergent_jacobi;
    conv_gauss_seidel(i) = convergent_gauss_seidel;

end

%{ 

    Part c)

        Plotting results

%}

% Iteration plot
figure;
plot(iter_jacobi, 'DisplayName', 'Jacobi');
hold on;
plot(iter_gauss_seidel, 'DisplayName', 'Gauss-Seidel');
% Adding tittle
title('Distribution of number of iterations');
% Adding legen
legend;
% Adding x-label
xlabel('Number of trials');
% Adding y-label
ylabel('Number of iterations');
hold off;

% Convergence plot
figure;
bar([sum(conv_jacobi), sum(conv_gauss_seidel)]);
set(gca, 'XTickLabel', {'Jacobi', 'Gauss-Seidel'});
title('Number of convergent trials');
ylabel('Number of trials');

%{ 

    Part d)

        Answering questions

%}

% Question & Answer 1
fprintf('\n')
fprintf('Which method seems more reliable? \n')
fprintf(['Answer: Given the tests done, in terms of reliability, both methods have similar requirements. They ' ...
    'both require the matrix to have certain properties (like diagonal dominance) ' ...
    'for guaranteed convergence. However, Gauss-Seidel might be considered slightly ' ...
    'more reliable in practice since it often converges faster and might do so in ' ...
    'situations where Jacobi might take much longer. \n'])
fprintf('\n')

% Question & Answer 2
fprintf('\n')
fprintf('Which method requires fewer iterations? \n')
fprintf(['Answer: Gauss-Seidel method requires fewer iterations than the ' ...
    'Jacobi method to converge to a solution, especially for diagonally dominant systems. In the Jacobi method, ' ...
    'all updates for the current iteration are calculated using the values from the previous iteration. ' ...
    'In contrast, the Gauss-Seidel method updates the solution values immediately as they are computed. This ' ...
    'means that Gauss-Seidel can utilize already updated values from the current iteration to compute ' ...
    'the subsequent ones, which can accelerate convergence. \n'])
fprintf('\n')

%{ 

    Part e)

        Repeating the experiment with AT*A=b

%}

% Initializing n an number of trials variables
n2 = 100;
num_trials2 = 1000;

% Initializing that will store the number of iterations each method took 
% for each trial and whether or not they converged.
iter_jacobi2 = zeros(num_trials2, 1);
iter_gauss_seidel2 = zeros(num_trials2, 1);
conv_jacobi2 = zeros(num_trials2, 1);
conv_gauss_seidel2 = zeros(num_trials2, 1);

% For-loop for trials
for i = 1:num_trials2

    % Generating a random matrix A
    A2 = randn(n, n);

    % While-loop checking the rank of A. If A is singular(rank less than n)
    % a new random matrix is generated until get a non-singular matrix(full
    % rank matrix).
    while rank(A2) < n

        A2 = randn(n2, n2);

    end

    % Generating a random vector b
    b2 = randn(n2, 1);

    % Calling the jacobi and gauss_seidel functions implementations to 
    % solve the linear system
    [~, niter_jacobi2, convergent_jacobi2] = jacobi(A2' * A2, A2' * b2, 1e-6);
    [~, niter_gauss_seidel2, convergent_gauss_seidel2] = gauss_seidel(A2' * A2, A2' * b2, 1e-6);

    % Storing the results of each trial in the arrays initialized 
    iter_jacobi2(i) = niter_jacobi2;
    iter_gauss_seidel2(i) = niter_gauss_seidel2;
    conv_jacobi2(i) = convergent_jacobi2;
    conv_gauss_seidel2(i) = convergent_gauss_seidel2;

end

%{ 

    Part e.1)

        Plotting results

%}

% Iteration plot
figure;
plot(iter_jacobi2, 'DisplayName', 'Jacobi');
hold on;
plot(iter_gauss_seidel2, 'DisplayName', 'Gauss-Seidel');
% Adding tittle
title('Distribution of number of iterations with AT*A=b.');
% Adding legend
legend;
% Adding x-label
xlabel('Number of trials');
% Adding y-label
ylabel('Number of iterations');
hold off;

% Convergence plot
figure;
bar([sum(conv_jacobi2), sum(conv_gauss_seidel2)]);
set(gca, 'XTickLabel', {'Jacobi', 'Gauss-Seidel'});
title('Number of convergent trials with AT*A=b.');
ylabel('Number of trials');

%{ 

    Part e.2)

        Answering questions

%}

% Question & Answer 1
fprintf('\n')
fprintf('Which method seems more reliable? \n')
fprintf(['Answer: Given that AT*A is symmetric and positive definite (provided  A is non-singular), ' ...
    'the Gauss-Seidel method is generally more reliable than the Jacobi method for convergence. This means ' ...
    'that for a symmetric and positive definite matrix, Gauss-Seidel will typically converge to the solution ' ...
    'more often than Jacobi. \n'])
fprintf('\n')

% Question & Answer 2
fprintf('\n')
fprintf('Which method requires fewer iterations? \n')
fprintf(['Answer: Gauss-Seidel method requires fewer iterations than the Jacobi method to converge ' ...
    'to a solution.  It is the same case that Ax=b, even with fewer iterations. \n'])
fprintf('\n')
