function test_suite=test_calc_frac_grad_expand_mat
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
    assertExceptionThrown(@()calc_frac_grad_expand_mat('foobar'), '*');
end

function test_frac_grad_expand

Nf_f = 19;
N_fractures = 4;
Nf_i = [6 5 3 5];

input = (1:Nf_f)';

%% Expanding gradient matrix
frac_grad_expand_mat = calc_frac_grad_expand_mat(N_fractures,Nf_f, Nf_i);

output = frac_grad_expand_mat*input;
assertEqual(output,[0;1;1;1;1;1;0;0;1;1;1;1;0;0;1;1;0;0;1;1;1;1;0]);

end


