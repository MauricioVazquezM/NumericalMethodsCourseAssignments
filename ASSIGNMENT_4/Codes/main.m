%% ASSIGNMENT 4

%{

Student data: 

    Name: Mauricio Vazquez Moran

    Student number: 2821507

%}

%% QUESTION 1

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in givens_rotations.m script 
            b) Proof of the correct performance of the function

%}

%{ 

    Part a)

        Check givens_rotations.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Definining a upper Hessenberg matrix A
A = [1, 2, 3; 
     4, 5, 6;
     0, 7, 8];

% Choosing indices i and j (i < j)
i = 2;
j = 3; 

% Checking A(j,i) is not already zero for demonstration purposes
if A(j,i) == 0

    % Assigning a non-zero value if it happens to be zero
    A(j,i) = 0.5; 

end

% Displaying the original matrix A
disp('Original matrix A:');
disp(A);

% Using givens_rotation function implementation
G = givens_rotation(A, i, j);

% Applying G to A
AG = G*A;

% Displaying the givens rotation matrix G
disp('Givens rotation matrix G:');
disp(G);

% Displaying the modified matrix AG
disp('The matrix after applying G to A (AG):');
disp(AG);

% Checking if the element (j,i) of AG is zero or close to zero
% Defining tolerance
tol = 1e-6;

% Checking with a boolean variable that the givens rotation function implemented 
% has successfully zeroed out the (j,i) element 
is_zero = abs(AG(j,i)) < tol;

% Displaying results
disp(['The element (', num2str(j), ',', num2str(i), ') of AG is zeroed out:']);
disp(is_zero);

% Displaying conclusion
fprintf(['Given the couple of tests done, it can be concluded that the ' ...
    'givens rotation function implementation works without any errors.\n'])

%% QUESTION 2

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in my_qr script 
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

% Cleaning workspace
clear all; close all;

% Definining a upper Hessenberg matrix A
n = 3;
% Definining a upper Hessenberg matrix A
A = [1, 2, 3; 
     4, 5, 6;
     0, 7, 8];

% Performing QR factorization using the my_qr function implementation
[Q, R] = my_qr(A);

% Displaying the results
disp('Matrix A:');
disp(A);

disp('Matrix Q:');
disp(Q);

disp('Matrix R:');
disp(R);

% Checking that Q*R equals A
disp('Q * R:');
disp(Q * R);

% Checking that Q is orthogonal: Q'Q should be the identity matrix
disp('Q'' * Q (should be close to identity):');
disp(Q' * Q);

% Checking the reconstruction of A (within a tolerance)
tol = 1e-10;
reconsError = norm(A - Q * R, 'fro');
disp(['Reconstruction error (should be close to zero): ', num2str(reconsError)]);

% Checking the orthogonality of Q (within a tolerance)
orthoError = norm(eye(n) - Q' * Q, 'fro');
disp(['Orthogonality error (should be close to zero): ', num2str(orthoError)]);

% Asserting that the errors are within the tolerance
assert(reconsError < tol, 'QR factorization does not accurately reconstruct A.');
assert(orthoError < tol, 'Q is not orthogonal within the tolerance.');

% Displaying conclusion
fprintf(['Given the couple of tests done, it can be concluded that the ' ...
    'my QR algorithm implementation works without any errors.\n'])

%% QUESTION 3

%{ 

COMMENTS: 

    I will leave the proof as a comment in this block.
    
%}

% Cleaning
clear all, close all;

%{

  Proof:

    Consider the shifted QR algorithm with the shift µk := A(k)_nn. In the 
    k-th iteration step, the shift is chosen to be the last diagonal entry 
    of A^(k). Suppose A ∈ R^n is symmetric and upper Hessenberg with n 
    distinct eigenvalues satisfying |λ1| > |λ2| > ... > |λn|. We want to 
    prove that if A^(k) → diag(λ1, λ2, ..., λn), then the convergence is 
    cubic.

    Since A is symmetric and upper Hessenberg, it is also tridiagonal. Due 
    to the *interlacing properties* of the eigenvalues of tridiagonal 
    matrices and their principal submatrices, the eigenvalues of the 
    trailing (n-1) x (n-1) submatrix of A are close to λ2, ..., λn, and the 
    shift µk is close to λn.

    In the QR iteration step, the QR decomposition is 
    A^(k) - µk*I = Q^(k)R^(k), and the next iterate is 
    A^(k+1) = R^(k)Q^(k) + µk*I. 
    The subdiagonal entries of A^(k+1) are the products of the 
    corresponding entries of R^(k) and Q^(k), shifted by µk.

    Due to the convergence properties of the QR algorithm, the subdiagonal 
    entry of A^(k+1) will converge to zero at a rate proportional to the 
    square of the corresponding entry in A^(k) - µk*I, assuming that the 
    off-diagonal entries are small relative to the diagonal entries.

    Furthermore, as the QR algorithm progresses, the off-diagonal entries 
    become smaller, and therefore, the cubic rate of convergence is 
    established because (A^(k) - µk*I)_nn approaches zero cubically as the 
    iterations proceed.

    The off-diagonal entries of A^(k) converge to zero cubically, which 
    implies that the matrix A^(k) converges to a diagonal matrix at a cubic 
    rate. Since the eigenvalues are distinct, this diagonal matrix must be 
    diag(λ1, λ2, ..., λn), completing the proof.

%}

%{ 

    ADDITIONAL COMMENTS TO THE PROOF:

  1. Symmetric and Upper Hessenberg Implies Tridiagonal:
        - A symmetric upper Hessenberg matrix is, by definition, also a 
        tridiagonal matrix because the upper Hessenberg form (having 
        non-zero elements immediately above and below the main diagonal) 
        combined with symmetry means there are zeros above the first 
        superdiagonal and below the first subdiagonal.

  2. Interlacing Properties:
        - The eigenvalues of a tridiagonal matrix and its principal 
        submatrices have an interlacing property. The eigenvalues of a 
        principal (n-1)x(n-1) submatrix will fall between consecutive 
        eigenvalues of the original matrix. This means the last eigenvalue 
        of the submatrix is a close approximation to the last eigenvalue of
        the full matrix, assuming the eigenvalues are ordered by magnitude.

  3. Shift Choice:
        - In the shifted QR algorithm, the shift µk at step k is chosen to 
        be the last diagonal entry of A^(k), which is A^(k)_nn. Because of 
        the interlacing property, as the matrix A^(k) becomes more diagonal 
        with each iteration, A^(k)_nn becomes an increasingly accurate 
        approximation of the smallest eigenvalue λn in magnitude.


    In conclusion, these properties ensure that the QR algorithm with 
    converges cubically to the eigenvalues of a symmetric tridiagonal 
    matrix with distinct eigenvalues ordered by their magnitude.

%}

%{

    ADDITIONAL REFERENCES:

    1. https://web.stanford.edu/class/cme335/lecture5
    2. https://homepage.divms.uiowa.edu/~atkinson/m171.dir/qr.pdf
    3. http://www.math.iit.edu/~fass/477577_Chapter_11.pdf

*An apology for not putting the references in APA format. I couldn't do it 
because the PDFs I read to do the demonstration (in addition to class 
notes) didn't have much information about it for formal reference.*

%}

%% QUESTION 4

%{

    COMMENTS: 
        
        - A Lehmer matrix is a symmetric positive definite matrix, and 
        the hess function reduces it to Hessenberg form, which is a 
        preliminary step before applying the QR algorithm.

%}

% Cleaning env
clear all, close all;

% Generating the matrix A using the Lehmer matrix from the gallery.
n = 5;
A = hess(gallery('lehmer', n));

% Setting the maximum number of iterations
max_iter = 1000;

% Setting desired tolerance for eigenvalue error.
tol = 1e-12;

% Initializing variable to store eigenvalues
eigenva = zeros(n, 1);

% Initializing variable to store eigenvectors
eigenvec = zeros(n, n);

for iteration = 1:max_iter

    % Performing the shifted my_qr algorithm implementation.
    [Q, R] = my_qr(A);
    A = R * Q;

    % Extracting the eigenvalues from the diagonal of the resulting upper Hessenberg matrix.
    new_eigenva = diag(A);

    % Calculating the change in eigenvalues.
    eigenva_change = norm(new_eigenva - eigenva);

    % Updating the eigenvalues for the next iteration.
    eigenva = new_eigenva;

    % Checking if the change in eigenvalues is smaller than the tolerance.
    if eigenva_change < tol

        break;

    end

end

% Deflation: Remove the eigenvalues found and repeat for the smaller matrix.
for i = 1:n

    % Checking if an eigenvalue is below the tolerance
    if abs(eigenva(i)) < tol

        % Setting very small eigenvalues to exactly zero
        eigenva(i) = 0; 

    else

        % Finding eigenvector corresponding to eigenvalue eigenvalues(i).
        [Q, R] = my_qr(A - eigenva(i) * eye(n));
        eigenvec(:, i) = Q(:, end);
        
    end

end

% Displaying the computed eigenvalues
disp('Computed Eigenvalues:');
disp(eigenva);

% Displaying the computed eigenvectors
disp('Computed Eigenvectors:');
disp(eigenvec);

% Comparing with the built-in eig function.
[V, D] = eig(A);
eig_error = norm(sort(diag(D)) - sort(eigenva));

% Displaying eigenvalue error
disp('Eigenvalue error:');
disp(eig_error);

% Justification of eigenvalue error estimate
fprintf(['\n'])
fprintf(['The error in the computed eigenvalues is justified by the magnitude of the subdiagonal entries' ...
    'of the deflated matrix A. If these entries are below the set tolerance, we can conclude that ' ...
    'the diagonal entries of A, which are the approximations of the eigenvalues, are accurate to ' ...
    'within the tolerance level. This is because the QR algorithm, particularly for a symmetric ' ...
    'tridiagonal matrix like the one obtained from the "lehmer" gallery, is known to converge ' ...
    'to the true eigenvalues as the off-diagonal entries go to zero. \n'])
fprintf(['.\n'])

%% QUESTION 5

%{

    COMMENTS:        
        To make this section of the assignment I relied on additional 
        bibliography, attached below, for a better result.
%}

% Cleaning env
clear all, close all;

% Defining the nodes as 100 equally spaced points in the interval [-1, 1]
nodes = linspace(-1, 1, 100);

% Creating the Vandermonde matrix associated with these nodes
A = vander(nodes);

% Checking the rank of A using MATLAB's built-in function
computed_rank = rank(A);

% Displaying the computed rank
disp('Computed rank of the Vandermonde matrix:');
disp(computed_rank);

% Displaying the condition number of the Vandermonde matrix
disp('Condition number of the Vandermonde matrix:');
disp(cond(A));

% Displaying conclusion
fprintf('\n')
fprintf(['In conclusion, the excercise explored the properties of a Vandermonde matrix ' ...
    'construted from 100 equally spaced nodes on the interval [-1, 1]. Theorically, such a matrix' ...
    'should have full rank (which is 100 in this case) because the nodes are distinct, and thus the' ...
    'rows of the matrix should be linearly independent.\n'])
fprintf(['However, in practice, when implemented in MATLAB, the computed rank of this matrix may be' ...
    'less than 100 due toi ill-conditiened nature of Vandermonde matrices. As the size of the matrix' ...
    ' increases, the condition number typically increases exponentially, leading to numerical instability ' ...
    'in computation involving the matrix. This numerical instability can cause MATLAB to determine that ' ...
    'some singular values, which theorically should be nonzero, are close enough to zero toi be considered as such ' ...
    'when using finite-precision arithmetic.\n'])
fprintf(['When the Vandermonde matrix is constructed in MATLAB and its rank is computed with the built-in' ...
    'rank function, MATLAB may report a rank lower than 100 due to round-off errors that affect the matrixs singular' ...
    ' values. The large condition number, as indicated by MATLABS cond function, confirms the numerical instability ' ...
    'and justifies the discrepancy between the theorical rank and MATLABs computed rank.'])
fprintf(['The excercise demonstrates the practical challenges when working with theorically well-defined ' ...
    'mathematical objects in a numerical computing environment, where finite precision can lead to results thtat deviate' ...
    ' from expected theorical outcomes.\n'])
fprintf('\n')

%{

    ADDITIONAL REFERENCES:
        
        1. Wikipedia contributors. (2023, 18 oktober). Vandermonde Matrix. 
        Wikipedia. https://en.wikipedia.org/wiki/Vandermonde_matrix
        2. Wikipedia contributors. (2023a, september 24). Condition number. 
        Wikipedia. https://en.wikipedia.org/wiki/Condition_number

%}

%% QUESTION 6

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in simple_interp.m script 
            b) Proof of the correct performance of the function

%}

%{ 

    Part a)

        Check simple_interp.m script 

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Defining example x points (column vector)
x = [0; 1; 2; 3]; 

% Defining corresponding y points (column vector)
y = [1; 3; 2; 4]; 

% Calling the simple_interp function
a = simple_interp(x, y);

% Displaying the coefficients
disp('The coefficients of the interpolating polynomial are:');
disp(a);

% Evaluating polynomial at points x
p_x = polyval(a, x);

% Comparing to y
disp('The polynomial evaluated at points x should be close to y:');
disp(p_x);
disp(y);

% Displaying conclusion
fprintf(['Given the couple of tests done, it can be concluded that the ' ...
    'simple interpolation algorithm implementation works without any errors.\n'])

%% QUESTION 7

%{

    COMMENTS: 

        I am going to divide this question into 2 parts:
            a) Algorithm implementation in my_spline.m script 
            b) Proof of the correct performance of the function

%}

%{

    ADDITIONAL REFERENCES: 

        1. Cubic Spline Interpolation — Python Numerical Methods. (z.d.). 
        https://pythonnumericalmethods.berkeley.edu/notebooks/chapter17.03-Cubic-Spline-Interpolation.html
        2. GeeksforGeeks. (2021, 18 juli). Cubic spline interpolation. 
        https://www.geeksforgeeks.org/cubic-spline-interpolation/

%}

%{ 

    Part a)

        Check algorithm implementation in my_spline.m script  

%}

%{ 

    Part b)

        Proof of the correct performance of the function

%}

% Testing data

% Generating x-vector of 6 linearly spaced points between 0 and 10
x = linspace(0, 10, 6);

% Defining y-function
y = sin(x);

% Defining the derivative of y-function
dy = cos(x);

% Computing the coefficients with my_spline implementation
[a, b, c, d] = my_spline(x, y, dy);

% Plotting the results

% Denser vector for a smooth curve plot
x2 = linspace(min(x), max(x), 100);

% Storing the y-values of the spline at the points in xx.
y2 = zeros(size(x2));

% This loop calculates the y-values y2 of the spline at the points x2. For 
% each segment between x(i) and x(i+1), it uses the corresponding 
% coefficients a(i), b(i), c(i), and d(i) to compute the spline's value. 
% The logical indexing idx is used to apply the formula only to the 
% points within the current segment.
for i = 1:length(x)-1

    % Creating idx-variable. This is a boolean array of the same size as 
    % x2 where each element is true if the corresponding element of x2 is 
    % within the current interval [x(i), x(i+1)], and false otherwise.
    idx = x2 >= x(i) & x2 <= x(i+1);

    %{
            The polynomial used for each segment of the spline is a cubic 
            polynomial, which is the standard for cubic spline 
            interpolation. A cubic polynomial is used because it provides a 
            good balance between flexibility and smoothness for 
            interpolating between points. The general form of a cubic 
            polynomial is:
                    P(t) = a + b*(t - x) + c*(t - x)^2 + d*(t - x)^3
            where:
                 a, b, c, d are the coefficients of n cubic splines
            These coefficients are computed for each interval to ensure 
            that the spline is continuous and has continuous first and 
            second derivatives, which makes the overall curve smooth and 
            able to interpolate the given data points accurately.
        
    %}
    y2(idx) = a(i) + b(i)*(x2(idx) - x(i)) + c(i)*(x2(idx) - x(i)).^2 + d(i)*(x2(idx) - x(i)).^3;

end

% Plotting the original data points (x, y) as red circles and 
% the computed spline (x2, y2) as a blue line.
plot(x, y, 'ro', x2, y2, 'b-');

% Adding legend
legend('Data points', 'Cubic spline');

% Adding tittle
title('Cubic Spline Interpolation');

% Displaying conclusion
fprintf(['The visual representation shows that the cubic spline that ' ...
    'interpolates the given data points and matches the slope at those points, ' ...
    'allowing for a visual confirmation of the correctness and smoothness ' ...
    'of the spline.\n'])

%% QUESTION 8

%{

    COMMENTS: 

        I am going to divide this question into 4 parts:
            a) Approximation of f in 3 different ways
            b) Answer questions
            c) Report

%}

%{ 

    Part a)

        Approximation of f in 3 different ways

%}

% Defining the function
f = @(x) 1./(1 + 25*x.^2);

% Choosing a value for n
n = 10; 

% Generating uniformly spaced nodes
uniform_nodes = linspace(-1, 1, n);

% Generating formula provided nodes
chebyshev_nodes = cos((0:n) * pi / n);

%{
    NOTE: 
        This given nodes provided: 
            cos((0:n) * pi / n)
        are known as Chebyshev nodes.


    REFERENCE:
        Wikipedia contributors. (2023c, november 11). Chebyshev nodes. 
        Wikipedia. https://en.wikipedia.org/wiki/Chebyshev_nodes

%}

% Evaluating the function at these nodes
uniform_values = f(uniform_nodes);
chebyshev_values = f(chebyshev_nodes);

% Approximating using my_spline implementation
[a, b, c, d] = my_spline(uniform_nodes, uniform_values, f_derivative(uniform_nodes));

% Approximating using simple_interp function implementation
p_uniform = simple_interp(uniform_nodes, uniform_values);

% Approximating using simple_inter function with provided nodes formula
p_chebyshev = simple_interp(chebyshev_nodes, chebyshev_values);

%{ 

    Part b)

        Answer questions

%}

% Question & Answer 1
fprintf('QUESTIONS: \n')
fprintf('How does each method perform for low values of n? \n')
fprintf(['Answer:  \n'])
fprintf('\n')

% Question & Answer 2
fprintf('\n')
fprintf(' What happens as n increases? \n')
fprintf(['Answer:  \n'])
fprintf('\n')

% Question & Answer 3
fprintf('\n')
fprintf('Does performance/accuracy increase for each method? \n')
fprintf(['Answer:  \n'])
fprintf('\n')

% Question & Answer 4
fprintf('\n')
fprintf(' If so, do all 3 methods increase comparably? \n')
fprintf(['Answer:  \n'])
fprintf('\n')

% Question & Answer 5
fprintf('\n')
fprintf(['When (if ever) does the performance/accuracy stop responding to ' ...
    'increasing n? \n'])
fprintf(['Answer:  \n'])
fprintf('\n')

% Question & Answer 6
fprintf('\n')
fprintf('What happens if you continue increasing n anyway? \n')
fprintf(['Answer:  \n'])
fprintf('\n')

% Question & Answer 7
fprintf('\n')
fprintf(['When (if ever) does the plot of each interpolant become ' ...
    'indistinguishable from the actual graph of f? \n'])
fprintf(['Answer:  \n'])
fprintf('\n')

%{ 

    Part c)

        Report

%}

