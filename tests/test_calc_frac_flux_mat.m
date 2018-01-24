function test_suite=test_calc_frac_flux_mat
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
    assertExceptionThrown(@()calc_frac_flux_mat('foobar'), '*');
end

function test_frac_fluxes

Nf_f = 19;
N_fractures = 4;
Nf_i = [6 5 3 5];

Ni=zeros(size(Nf_i));
for i=1:N_fractures
    for j=i:-1:1
        Ni(i) = Ni(i) +Nf_i(j);
    end
end
Nii = 1:1:length(Ni);
Nf_ff = Nf_f +N_fractures;

[frac_left_flux_mat, frac_right_flux_mat] = calc_frac_flux_mat(N_fractures,Nf_f, Nf_i);

input = (1:Nf_ff)';

%% Test left fluxes
grad_left = frac_left_flux_mat*input;
assertEqual(grad_left,[0;2;3;4;5;6;0;9;10;11;12;0;15;16;0;19;20;21;22])

%% Test right fluxes
grad_right = frac_right_flux_mat*input;
assertEqual(grad_right,[2;3;4;5;6;0;9;10;11;12;0;15;16;0;19;20;21;22;0]);

end

function test_empty()
    % Test empty vector input
    assertEqual(calc_frac_flux_mat([],[],[]), sparse([]));
end


