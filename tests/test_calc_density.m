function test_suite=test_calc_density
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
    assertExceptionThrown(@()calc_density('foobar'), '*');
end

function test_scalar
[density_m, density_f]  = calc_density(1,1,1,1);
assertElementsAlmostEqual([round(density_m *1e4)/1e4 round(density_f *1e4)/1e4],[999.3311 999.3311]);
end

function test_vector
[density_m, density_f] = calc_density(ones(10,10),ones(10,1),ones(10,10),ones(10,1));
assertEqual([round(density_m *1e4)/1e4 round(density_f *1e4)/1e4],[999.3311*ones(10,10) 999.3311*ones(10,1)]);
end

function test_empty()
    % Test empty vector input
    assertEqual(calc_density([],[],[],[]), [[] [] [] []]);
end