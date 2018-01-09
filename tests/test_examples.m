function test_suite=test_examples
try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions=localfunctions();
catch % no problem; early Matlab versions can use initTestSuite fine
end
initTestSuite;
end

function test_ex1
clear_workspace()
load comsol_bench_h.mat

global k_ratio
k_ratio = 1e3; % Run with Kf/Km = 1e3

THERMAID('Input_ex1',0)

udata = evalin('base','udata');
x     = evalin('base','x');
p     = evalin('base','p');
pf    = evalin('base','pf');

xf = 0:udata.dxf:5;
xf = xf(1:udata.Nf_i(1));

pm_1e3_EDFM   = p(:,floor(udata.Nf(2)/2));
pf_1e3_EDFM   = pf(udata.Nf_i(1)+1:udata.Nf_f);

pf_interp_1e3 = interp1(x_pf_1e3,pf_1e3,xf');
pm_interp_1e3 = interp1(x_pm_1e3,pm_1e3,x');
pm_interp_1e3(isnan(pm_interp_1e3)) = 0;

RMSE_pf_1e3 = sqrt(mean((pf_1e3_EDFM-pf_interp_1e3).^2))/(max(pf_interp_1e3)-min(pf_interp_1e3));
RMSE_pm_1e3 = sqrt(mean((pm_1e3_EDFM-pm_interp_1e3).^2))/(max(pm_interp_1e3)-min(pm_interp_1e3));

assertEqual(round(RMSE_pf_1e3 *1e7)/1e7,0.0050262);
assertEqual(round(RMSE_pm_1e3 *1e7)/1e7,0.0046731);

end

function test_ex2
clear_workspace()
load comsol_bench_th2.mat

global k_ratio

% Run with Kf/Km = 1e5
k_ratio = 1e5;
THERMAID('Input_ex2',0)

udata = evalin('base','udata');
XY1     = evalin('base','XY1');
tNewf    = evalin('base','tNewf');

%% Post analysis
xf = udata.dxf/2:udata.dxf:abs(XY1(1,2)-XY1(end,3));

T_horz = tNewf(udata.Nf_i(1)+1:udata.Nf_f);
T_vert = tNewf(1:udata.Nf_i(1));

%% Quantitative analysis
x_interp = xf';
T_interp_horz = interp1(x_ref_T_horz,T_ref_horz,x_interp);

% Remove NaN from interpolation
T_interp_horz(end) = T_interp_horz(end-2);
T_interp_horz(end-1) = T_interp_horz(end-2); 

T_interp_vert = interp1(x_ref_T_vert,T_ref_vert,x_interp);
T_interp_horz(end) = T_interp_horz(end-1);
T_interp_vert(end) = T_interp_vert(end-1);

RMSE_T_vert = sqrt(mean((T_vert-T_interp_vert).^2))/(max(T_interp_vert)-min(T_interp_vert));
RMSE_T_horz = sqrt(mean((T_horz-T_interp_horz).^2))/(max(T_interp_horz)-min(T_interp_horz));

assertEqual(round(RMSE_T_vert *1e6)/1e6,0.052653);
assertEqual(round(RMSE_T_horz *1e6)/1e6,0.011209);
end

function test_ex3
clear_workspace()
global k_ratio
k_ratio = 1e5; % Run with Kf/Km = 1e5

THERMAID('Input_ex3',0);

pf    = evalin('base','pf');
p    = evalin('base','p');
tNewf    = evalin('base','tNewf');
tNew    = evalin('base','tNew');

load reference_ex3.mat

assertElementsAlmostEqual(p,p_ref)
assertElementsAlmostEqual(pf,pf_ref)
assertElementsAlmostEqual(tNew,T_ref)
assertElementsAlmostEqual(tNewf,Tf_ref)
end

function test_ex4
clear_workspace()

THERMAID('Input_ex4',0);

pf    = evalin('base','pf');
p    = evalin('base','p');

load reference_ex4.mat

assertElementsAlmostEqual(p,p_ref)
assertElementsAlmostEqual(pf,pf_ref)
end

function test_ex5
clear_workspace()

THERMAID('Input_ex5',0);

pf    = evalin('base','pf');
p    = evalin('base','p');
tNewf    = evalin('base','tNewf');
tNew    = evalin('base','tNew');

load reference_ex5.mat

assertElementsAlmostEqual(p,p_ref)
assertElementsAlmostEqual(pf,pf_ref)
assertElementsAlmostEqual(tNew,T_ref)
assertElementsAlmostEqual(tNewf,Tf_ref)

end

function clear_workspace()
clear all
end
