function test_suite=test_calc_frac_mean_mat
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
    assertExceptionThrown(@()calc_frac_mean_mat('foobar'), '*');
end

function test_frac_mean

Nf_f = 19;
N_fractures = 4;
Nf_i = [6 5 3 5];

input = 0.5*ones(Nf_f,1);
input(1:6) = 3;
input(7:11) = 6;
input(12:14) = 9;

frac_mean_mat = calc_frac_mean_mat(N_fractures,Nf_f, Nf_i);

%% Averaging Matrix
arit_mean = 0.5*frac_mean_mat*input;
assertEqual(arit_mean,[3;3;3;3;3;3;3;6;6;6;6;6;6;9;9;9;9;0.5;0.5;0.5;0.5;0.5;0.5]);

%% Harmonic Averaging Matrix
harm_mean = 2./(frac_mean_mat*(1./input));
assertEqual(harm_mean,[3;3;3;3;3;3;3;6;6;6;6;6;6;9;9;9;9;0.5;0.5;0.5;0.5;0.5;0.5]);

end



