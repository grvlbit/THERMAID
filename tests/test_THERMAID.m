function test_suite=test_THERMAID
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_suite_name()
name='THERMAID test suite';
suite=MOxUnitTestSuite(name);

assertEqual(name,getName(suite));
end