function test_suite=test_calc_d_mean
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
    assertExceptionThrown(@()calc_d_mean('foobar'), '*');
end

function test_scalar
    assertEqual(round(calc_d_mean(1,1)*1e4)/1e4 ,0.2357);
end

function test_empty()
    % Test empty vector input
    assertExceptionThrown(@()calc_d_mean([],[]), '*');
end
