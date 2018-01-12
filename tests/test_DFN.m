function test_suite=test_DFN
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_frac_complex_n13
udata.len     = [10 10];              % physical length of the domain in x and y direction [m]
udata.Nf      = [100 100];  	      % number of cells in x and y direction
udata.dx      = udata.len./udata.Nf;  % cell length [m]

frac_complex_n13

assertTrue(udata.dxf > min(udata.dx))
assertEqual(udata.N_fractures,13)
assertEqual(udata.Nf_f,390)
assertEqual(length(udata.frac_angle),udata.Nf_f)
assertEqual(size(XY1,1),udata.Nf_f)
assertEqual(udata.Nf_i,[70 18 49 28 42 31 21 24 22 28 8 28 21])
end

function test_frac_complex_n37
udata.len     = [10 10];              % physical length of the domain in x and y direction [m]
udata.Nf      = [100 100];  	      % number of cells in x and y direction
udata.dx      = udata.len./udata.Nf;  % cell length [m]

frac_complex_n37

assertTrue(udata.dxf > min(udata.dx))
assertEqual(udata.N_fractures,37)
assertEqual(udata.Nf_f,572)
assertEqual(length(udata.frac_angle),udata.Nf_f)
assertEqual(size(XY1,1),udata.Nf_f)
assertEqual(size(udata.Nf_i,2),udata.N_fractures)
end

function test_frac_cross
udata.len     = [10 10];              % physical length of the domain in x and y direction [m]
udata.Nf      = [100 100];  	      % number of cells in x and y direction
udata.dx      = udata.len./udata.Nf;  % cell length [m]

frac_cross

assertTrue(udata.dxf > min(udata.dx))
assertEqual(udata.N_fractures,2)
assertEqual(udata.Nf_f,100)
assertEqual(length(udata.frac_angle),udata.Nf_f)
assertEqual(size(XY1,1),udata.Nf_f)
assertEqual(size(udata.Nf_i,2),udata.N_fractures)
end

function test_frac_cross_rot45
udata.len     = [100 100];              % physical length of the domain in x and y direction [m]
udata.Nf      = [100 100];  	      % number of cells in x and y direction
udata.dx      = udata.len./udata.Nf;  % cell length [m]

frac_cross_rot45

assertTrue(udata.dxf > min(udata.dx))
assertEqual(udata.N_fractures,2)
assertEqual(udata.Nf_f,88)
assertEqual(length(udata.frac_angle),udata.Nf_f)
assertEqual(size(XY1,1),udata.Nf_f)
assertEqual(size(udata.Nf_i,2),udata.N_fractures)
end

function test_frac_load_fracman

frac_load_fracman

assertTrue(udata.dxf > min(udata.dx))
assertEqual(udata.N_fractures,523)
assertEqual(udata.Nf_f,4846)
assertEqual(length(udata.frac_angle),udata.Nf_f)
assertEqual(size(XY1,1),udata.Nf_f)
assertEqual(size(udata.Nf_i,2),udata.N_fractures)
end

function test_frac_load_fracsim

frac_load_fracsim

assertTrue(udata.dxf > min(udata.dx))
assertEqual(udata.N_fractures,245)
assertEqual(udata.Nf_f,5315)
assertEqual(length(udata.frac_angle),udata.Nf_f)
assertEqual(size(XY1,1),udata.Nf_f)
assertEqual(size(udata.Nf_i,2),udata.N_fractures)
end

function test_frac_single
udata.len     = [10 10];              % physical length of the domain in x and y direction [m]
udata.Nf      = [100 100];  	      % number of cells in x and y direction
udata.dx      = udata.len./udata.Nf;  % cell length [m]

frac_single

assertTrue(udata.dxf > min(udata.dx))
assertEqual(udata.N_fractures,1)
assertEqual(udata.Nf_f,35)
assertEqual(length(udata.frac_angle),udata.Nf_f)
assertEqual(size(XY1,1),udata.Nf_f)
assertEqual(size(udata.Nf_i,2),udata.N_fractures)
end

function test_read_dfn_data_from_file
    assertExceptionThrown(@()read_dfn_data_from_file('foobar'),'*');
    assertExceptionThrown(@()read_dfn_data_from_file('foobar','fracMan',1,2,3,4),'*');
    
    filename = 'fracmanTest.csv';
    [x1,y1,x2,y2] = read_dfn_data_from_file(filename,'fracman',2,3,5,6);
    
    assertEqual(length(x1),6774)
    assertEqual(length(x2),6774)
    assertEqual(length(y1),6774)
    assertEqual(length(y2),6774)
    
    filename = 'fracSimTest.csv';
    [x1,y1,x2,y2] = read_dfn_data_from_file(filename,'fracsim',2,3,4,5);
    
    assertEqual(length(x1),513)
    assertEqual(length(x2),513)
    assertEqual(length(y1),513)
    assertEqual(length(y2),513)

end