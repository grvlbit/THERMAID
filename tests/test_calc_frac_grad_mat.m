function test_suite=test_calc_frac_grad_mat
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_input_handling()  %#ok<*DEFNU>
    % Check that an error is thrown for a bad input. The second
    % parameter is the expected error message, which in our case is
    % omitted because rectifiedCubic throws the error with a call to
    % `assert()`, not `error()`.
    % Throw an error if input is a string.
    assertExceptionThrown(@()calc_frac_grad_mat('foobar'), '*');
end

function test_frac_grad

Nf_f = 19;
N_fractures = 4;
Nf_i = [6 5 3 5];
Nf_ff = Nf_f +N_fractures;

input = ones(Nf_ff,1);
input(1:7) = 3;
input(8:13) = 6;
input(14:17) = 9;

%% Gradient Matrix
frac_grad_mat = calc_frac_grad_mat(N_fractures,Nf_f, Nf_i);

grad = frac_grad_mat*input;
assertEqual(grad,zeros(Nf_f,1));

end

function test_empty()
    % Test empty vector input
    assertEqual(calc_frac_grad_mat([],[],[]), sparse([]));
end


