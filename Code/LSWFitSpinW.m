% Fit to CAMEA data for exchange constants.

%% Global Fit: Develop SpinW model for Hamiltonian

clear, close all

% J1 = 2*-0.028; J2 = 2*0.152; J3 = 2*0.004; Dc = -0.0; % Nikotin 1969
% J1 = 2*-0.028; J2 = 2*0.152; J3 = 2*0.000; Dc = -0.091/(5/2)^2; % Okazaki 1963.
% J1 = 2*-0.0354; J2 = 2*0.1499; J3 = 2*0.000; Dc = -0.131/(5/2)^2; % CAMEA 2023.
J1 = -0.067718; J2 = 0.3022; J3 = -0.0043851; Dc = -0.026705; % Previous fit

% Instrumental E-resolution extracted from MJOLNIR for HHL along (001).
% Will be fit to polynomial in sw_instrument.
sigE = [0.3090, 0.0947; 1.0043, mean([0.0917, 0.0948]);...
    2.0002, mean([0.0889, 0.0919]); 2.9962, mean([0.1340, 0.1372]);...
    3.9921, mean([0.1321, 0.1351]); 5.0069, mean([0.1793, 0.1827]);...
    6.0028, mean([0.1779, 0.1811]); 6.9988, mean([0.2276, 0.2310]);...
    7.9947, mean([0.2264, 0.2298])];
fwhmE = sigE; fwhmE(:,2) = sigE(:,2)*2*sqrt(2*log(2));

MnF2Fit = spinw;
MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
    [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label',...
    {'MMn2' 'F'}, 'color', {'red' 'yellow'});

MnF2Fit.gencoupling('maxDistance', 10);
% MnF2Fit.table('bond',1:10)

MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
MnF2Fit.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))

MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
MnF2Fit.addaniso('Dc');

MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);

% Relax magnetic structure
% MnF2Fit.optmagsteep('nRun', 1e4)
% MnF2Fit.table('mag')

% plot(MnF2Fit, 'range', [1 1 1], 'bondMode', 'line', 'bondLinewidth0', 3)

% Global Fit: Prepare fitting to the dispersion (tutorial 35)

x1 = [J1 J2 J3 Dc];
% x1 = [J1 J2 Dc];
MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2' 'J3' 'Dc(3,3)'}, 'init', true)
% MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2' 'Dc(3,3)'}, 'init', true)

MnF2Fit.table('matrix')

pref = swpref;
pref.tid = 0;

par_fit          = struct;
par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';

par_fit.Evect    = linspace(0, 10, 1001);
par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                    'mat', {'J1' 'J2' 'J3' 'Dc(3,3)'}, 'init', true);
par_fit.xmin      = [-0.6 -0.6 -0.2 -0.5];
par_fit.xmax      = [0.6 0.6 0.2 0.0];
% par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
%                     'mat', {'J1' 'J2' 'Dc(3,3)'}, 'init', true);
% par_fit.xmin      = [-0.6 -0.6 -0.5];
% par_fit.xmax      = [0.6 0.6 0.0];
par_fit.x0        = x1;
par_fit.plot      = true;
par_fit.hermit    = false;
par_fit.optimizer = 'pso';
par_fit.maxiter   = 5e2;
par_fit.nrun      = 1;
par_fit.TolX      = 10^(-4);
par_fit.TolFun    = 10^(-5);
par_fit.iFact     = 0;

% Global Fit 1: Fit to the dispersion (tutorial 35)

fitStr = MnF2Fit.fitspec(par_fit);

disp(['p: ' num2str(fitStr.x)])
disp(['chi2r: ' num2str(fitStr.redX2)])

%% Global Fit 1: Extract uncertainties

% J1
delJ1 = -0.0024:0.0002:0.0024;
J1chi2 = zeros(size(delJ1));
J1Iter = zeros(size(delJ1));
T = sw_readtable('U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt');
dof = length(T) - 3 + 1; % https://spinw.org/spinwdoc/spinw_fitspec
for i=1:length(delJ1)
    J1 = fitStr.x(1) + delJ1(i); J2 = fitStr.x(2);...
        J3 = fitStr.x(3); Dc = fitStr.x(4);
    
    MnF2Fit = spinw;
    MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
        [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
    MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0],...
        'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
    MnF2Fit.gencoupling('maxDistance', 10);
    MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
    MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
    MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
    MnF2Fit.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
    MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
    MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
    MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
    MnF2Fit.addaniso('Dc', 1);
    MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);
    
    x1 = [J2 J3 Dc];
    MnF2Fit.matparser('param', x1, 'mat', {'J2' 'J3' 'Dc(3,3)'}, 'init', true)
    
    MnF2Fit.table('matrix')
    
    pref = swpref;
    pref.tid = 0;
    
    par_fit          = struct;
    par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';
    
    par_fit.Evect    = linspace(0, 10, 1001);
    par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                        'mat', {'J2' 'J3' 'Dc(3,3)'}, 'init', true);
    par_fit.xmin      = [-0.6 -0.2 -0.5];
    par_fit.xmax      = [0.6 0.2 0.0];
    par_fit.x0        = x1;
    par_fit.plot      = false;
    par_fit.hermit    = false;
    par_fit.optimizer = 'pso';
    par_fit.maxiter   = 2e2;
    par_fit.nrun      = 1;
    par_fit.TolX      = 10^(-4);
    par_fit.TolFun    = 10^(-5);
    
    fitStrIter = MnF2Fit.fitspec(par_fit);

    J1Iter(i) = fitStr.x(1) + delJ1(i);
    J1chi2(i) = fitStrIter.redX2;
    
    disp(['i: ' num2str(i)])
    disp(['p: ' num2str(fitStrIter.x)])
    disp(['chi2R: ' num2str(fitStrIter.redX2)])
end
chi2RThr = (1 + 1/dof)*min(J1chi2); % Threshold reduced chi-squared
indClose = find(J1chi2<chi2RThr); % J1's within the threshold
sigJ1 = (J1Iter(indClose(end)+1) - J1Iter(indClose(1)-1))/2; % Go to points just beyond intercept with threshold, subtract, divide by two to slightly overestimate uncertainty
figure
hold on
scatter(J1Iter, J1chi2)
yline(chi2RThr)
xlabel('J1 (meV)')
ylabel('Chi2R')
title(['J1: ' num2str(fitStr.x(1)) '. sigJ1: ' num2str(sigJ1)])
box on
hold off
saveas(gcf, strcat('U:/Data/Figures/MnF2/Fitting/UncJ1.pdf'))
disp(['J1: ' num2str(fitStr.x(1))])
disp(['sigJ1: ' num2str(sigJ1)])

% J2
delJ2 = -0.002:0.0001:0.002;
J2chi2 = zeros(size(delJ2));
J2Iter = zeros(size(delJ2));
for i=1:length(delJ2)
    J1 = fitStr.x(1); J2 = fitStr.x(2) + delJ2(i);...
        J3 = fitStr.x(3); Dc = fitStr.x(4);
    
    MnF2Fit = spinw;
    MnF2Fit.genlattice('lat_const', [4.87380 4.87380 3.31070], 'angled',...
        [90 90 90], 'sym', 136);
    MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0],...
        'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
    MnF2Fit.gencoupling('maxDistance', 10);
    MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
    MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
    MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
    MnF2Fit.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
    MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
    MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
    MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
    MnF2Fit.addaniso('Dc', 1);
    MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);
    
    x1 = [J1 J3 Dc];
    MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J3' 'Dc(3,3)'}, 'init', true)
    
    MnF2Fit.table('matrix')
    
    pref = swpref;
    pref.tid = 0;
    
    par_fit          = struct;
    par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';
    
    par_fit.Evect    = linspace(0, 10, 1001);
    par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                        'mat', {'J1' 'J3' 'Dc(3,3)'}, 'init', true);
    par_fit.xmin      = [-0.6 -0.2 -0.5];
    par_fit.xmax      = [0.6 0.2 0.0];
    par_fit.x0        = x1;
    par_fit.plot      = false;
    par_fit.hermit    = false;
    par_fit.optimizer = 'pso';
    par_fit.maxiter   = 2e2;
    par_fit.nrun      = 1;
    par_fit.TolX      = 10^(-4);
    par_fit.TolFun    = 10^(-5);
    
    fitStrIter = MnF2Fit.fitspec(par_fit);

    J2Iter(i) = fitStr.x(2) + delJ2(i);
    J2chi2(i) = fitStrIter.redX2;
    
    disp(['i: ' num2str(i)])
    disp(['p: ' num2str(fitStrIter.x)])
    disp(['chi2R: ' num2str(fitStrIter.redX2)])
end
chi2RThr = (1 + 1/dof)*min(J2chi2); % Threshold reduced chi-squared
indClose = find(J2chi2<chi2RThr); % J2's within the threshold
sigJ2 = (J2Iter(indClose(end)+1) - J2Iter(indClose(1)-1))/2; % Go to points just beyond intercept with threshold, subtract, divide by two to slightly overestimate uncertainty
figure
hold on
scatter(J2Iter, J2chi2)
yline(chi2RThr)
xlabel('J2 (meV)')
ylabel('Chi2R')
title(['J2: ' num2str(fitStr.x(2)) '. sigJ2: ' num2str(sigJ2)])
box on
hold off
saveas(gcf, strcat('U:/Data/Figures/MnF2/Fitting/UncJ2.pdf'))
disp(['J2: ' num2str(fitStr.x(2))])
disp(['sigJ2: ' num2str(sigJ2)])

% J3
delJ3 = -0.0015:0.0001:0.0015;
J3chi2 = zeros(size(delJ3));
J3Iter = zeros(size(delJ3));
for i=1:length(delJ3)
    J1 = fitStr.x(1); J2 = fitStr.x(2);...
        J3 = fitStr.x(3) + delJ3(i); Dc = fitStr.x(4);
    
    MnF2Fit = spinw;
    MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
        [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
    MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0],...
        'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
    MnF2Fit.gencoupling('maxDistance', 10);
    MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
    MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
    MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
    MnF2Fit.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
    MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
    MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
    MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
    MnF2Fit.addaniso('Dc', 1);
    MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);
    
    x1 = [J1 J2 Dc];
    MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2' 'Dc(3,3)'}, 'init', true)
    
    MnF2Fit.table('matrix')
    
    pref = swpref;
    pref.tid = 0;
    
    par_fit          = struct;
    par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';
    
    par_fit.Evect    = linspace(0, 10, 1001);
    par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                        'mat', {'J1' 'J2' 'Dc(3,3)'}, 'init', true);
    par_fit.xmin      = [-0.6 -0.6 -0.5];
    par_fit.xmax      = [0.6 0.6 0.0];
    par_fit.x0        = x1;
    par_fit.plot      = false;
    par_fit.hermit    = false;
    par_fit.optimizer = 'pso';
    par_fit.maxiter   = 2e2;
    par_fit.nrun      = 1;
    par_fit.TolX      = 10^(-4);
    par_fit.TolFun    = 10^(-5);
    
    fitStrIter = MnF2Fit.fitspec(par_fit);

    J3Iter(i) = fitStr.x(3) + delJ3(i);
    J3chi2(i) = fitStrIter.redX2;
    
    disp(['i: ' num2str(i)])
    disp(['p: ' num2str(fitStrIter.x)])
    disp(['chi2R: ' num2str(fitStrIter.redX2)])
end
chi2RThr = (1 + 1/dof)*min(J3chi2); % Threshold reduced chi-squared
indClose = find(J3chi2<chi2RThr); % J3's within the threshold
sigJ3 = (J3Iter(indClose(end)+1) - J3Iter(indClose(1)-1))/2; % Go to points just beyond intercept with threshold, subtract, divide by two to slightly overestimate uncertainty
figure
hold on
scatter(J3Iter, J3chi2)
yline(chi2RThr)
xlabel('J3 (meV)')
ylabel('Chi2R')
title(['J3: ' num2str(fitStr.x(3)) '. sigJ3: ' num2str(sigJ3)])
box on
hold off
saveas(gcf, strcat('U:/Data/Figures/MnF2/Fitting/UncJ3.pdf'))
disp(['J3: ' num2str(fitStr.x(3))])
disp(['sigJ3: ' num2str(sigJ3)])

% Dc
delDc = -0.003:0.0002:0.003;
Dcchi2 = zeros(size(delDc));
DcIter = zeros(size(delDc));
for i=1:length(delDc)
    J1 = fitStr.x(1); J2 = fitStr.x(2);...
        J3 = fitStr.x(3); Dc = fitStr.x(4) + delDc(i);
    
    MnF2Fit = spinw;
    MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
        [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
    MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0],...
        'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
    MnF2Fit.gencoupling('maxDistance', 10);
    MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
    MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
    MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
    MnF2Fit.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
    MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
    MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
    MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
    MnF2Fit.addaniso('Dc', 1);
    MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);
    
    x1 = [J1 J2 J3];
    MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2' 'J3'}, 'init', true)
    
    MnF2Fit.table('matrix')
    
    pref = swpref;
    pref.tid = 0;
    
    par_fit          = struct;
    par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';
    
    par_fit.Evect    = linspace(0, 10, 1001);
    par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                        'mat', {'J1' 'J2' 'J3'}, 'init', true);
    par_fit.xmin      = [-0.6 -0.6 -0.2];
    par_fit.xmax      = [0.6 0.6 0.2];
    par_fit.x0        = x1;
    par_fit.plot      = false;
    par_fit.hermit    = false;
    par_fit.optimizer = 'pso';
    par_fit.maxiter   = 2e2;
    par_fit.nrun      = 1;
    par_fit.TolX      = 10^(-4);
    par_fit.TolFun    = 10^(-5);
    
    fitStrIter = MnF2Fit.fitspec(par_fit);

    DcIter(i) = fitStr.x(4) + delDc(i);
    Dcchi2(i) = fitStrIter.redX2;
    
    disp(['i: ' num2str(i)])
    disp(['p: ' num2str(fitStrIter.x)])
    disp(['chi2R: ' num2str(fitStrIter.redX2)])
end
chi2RThr = (1 + 1/dof)*min(Dcchi2); % Threshold reduced chi-squared
indClose = find(Dcchi2<chi2RThr); % Dc's within the threshold
sigDc = (DcIter(indClose(end)+1) - DcIter(indClose(1)-1))/2; % Go to points just beyond intercept with threshold, subtract, divide by two to slightly overestimate uncertainty
figure
hold on
scatter(DcIter, Dcchi2)
yline(chi2RThr)
xlabel('Dc (meV)')
ylabel('Chi2R')
title(['Dc: ' num2str(fitStr.x(4)) '. sigDc: ' num2str(sigDc)])
box on
hold off
saveas(gcf, strcat('U:/Data/Figures/MnF2/Fitting/UncDc.pdf'))
disp(['Dc: ' num2str(fitStr.x(4))])
disp(['sigDc: ' num2str(sigDc)])

%% Global Fit: Construct the best-fit Hamiltonian

J1 = fitStr.x(1); J2 = fitStr.x(2); J3 = fitStr.x(3); Dc = fitStr.x(4);

MnF2Fit = spinw;
MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
    [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label',...
    {'MMn2' 'F'}, 'color', {'red' 'yellow'});
MnF2Fit.gencoupling('maxDistance', 10);
MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
MnF2Fit.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
MnF2Fit.addaniso('Dc', 1);
MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);

%% Global Fit: Read in the data

% Plot over the experimental data after loading the data.
T = sw_readtable('U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt');

% List of Q values from data.
Q = [[T(:).QH]; [T(:).QK]; [T(:).QL]]';
% List of energy values from data.
w = [T(:).EN1];
werr = [T(:).s1];
% List of intensity values from data.
inten = [T(:).I1];
% Indices corresponding to each slice
% HHL
ind1 = 1:5;
ind2 = 6:26;
ind3 = 27:40;
ind4 = 41:56;
ind5 = 57:67;
% H0L
ind6 = 68:84;
ind7 = 85:101;
ind8 = 102:112;
ind9 = 113:122;
% HK0
ind10 = 123:142;

save('U:\Data\MATLAB\MnF2\HamiltonianFit.mat')

%% Figure 1: Develop MnF2 Hamiltonian and simulate spectrum using SpinW

load('U:\Data\MATLAB\MnF2\HamiltonianFit.mat');
J1 = fitStr.x(1); J2 = fitStr.x(2); J3 = fitStr.x(3); J7 = 0.02; J8 = -0.02; Dc = fitStr.x(4); fwhmESplit = 0.1;

MnF2Splitting = spinw;
MnF2Splitting.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled', [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Splitting.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
MnF2Splitting.gencoupling('maxDistance', 10);
MnF2Splitting.addmatrix('label', 'J1', 'value', J1);
MnF2Splitting.addmatrix('label', 'J2', 'value', J2);
MnF2Splitting.addmatrix('label', 'J3', 'value', J3);
MnF2Splitting.addmatrix('label', 'J7', 'value', J7);
MnF2Splitting.addmatrix('label', 'J8', 'value', J8);
MnF2Splitting.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
MnF2Splitting.addcoupling('mat', 'J1', 'bond', 1);
MnF2Splitting.addcoupling('mat', 'J2', 'bond', 2);
MnF2Splitting.addcoupling('mat', 'J3', 'bond', 3);
MnF2Splitting.addcoupling('mat', 'J7', 'bond', 7);
MnF2Splitting.addcoupling('mat', 'J8', 'bond', 8);
MnF2Splitting.addaniso('Dc', 1);
MnF2Splitting.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);

% MnF2Splitting.energy
% MnF2Splitting.table('mag')
% plot(MnF2Splitting, 'range', [1 1 1], 'bondMode', 'line', 'bondLinewidth0', 3)
% MnF2Splitting.getmatrix('bond', 2);

spec = MnF2Splitting.spinwave({[-0.5,0.5,1] [0,0,1] [0.5,0.5,1] 223}, 'formfact', true, 'hermit', true);
spec = sw_neutron(spec, 'pol', true, 'n', [0 0 1]);
spec = sw_egrid(spec, 'Evect', linspace(0, 10, 1000), 'component', 'Pxy+Pxz', 'T', 0, 'ImagChk', true);

f = figure('Position', [2, 3.3, 15, 13.24], 'Units', 'centimeters');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile()
hold on
box on
axis square
xline(0)
yline(0)
plot([-0.5, 0, 0.5], [0.5, 0, 0.5], '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2);
text(-0.02, -0.1, '\Gamma', 'FontSize', 18, 'HorizontalAlignment', 'right');
text(0.55, 0.5, 'M', 'FontSize', 18);
text(0.55, -0.5, 'M', 'FontSize', 18);
text(-0.65, 0.5, 'M', 'FontSize', 18);
text(-0.65, -0.5, 'M', 'FontSize', 18);
text(0.55, 0, 'X', 'FontSize', 18);
text(0, 0.6, 'X', 'FontSize', 18, 'HorizontalAlignment', 'center');
text(-0.65, 0, 'X', 'FontSize', 18);
text(0, -0.6, 'X', 'FontSize', 18, 'HorizontalAlignment', 'center');
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 18)
xlim([-0.5, 0.5])
ylim([-0.5, 0.5])
hold off

nexttile()
hold on
box on
axis square
xline(0)
yline(0)
plot([-0.5, 0, 0.5], [0.5, 0, 0.5], '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2);
text(-0.02, -0.1, 'Z', 'FontSize', 18, 'HorizontalAlignment', 'right');
text(0.55, 0.5, 'A', 'FontSize', 18);
text(0.55, -0.5, 'A', 'FontSize', 18);
text(-0.65, 0.5, 'A', 'FontSize', 18);
text(-0.65, -0.5, 'A', 'FontSize', 18);
text(0.55, 0, 'R', 'FontSize', 18);
text(0, 0.6, 'R', 'FontSize', 18, 'HorizontalAlignment', 'center');
text(-0.65, 0, 'R', 'FontSize', 18);
text(0, -0.6, 'R', 'FontSize', 18, 'HorizontalAlignment', 'center');
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 18)
xlim([-0.5, 0.5])
ylim([-0.5, 0.5])
hold off

nexttile()
sw_plotspec(spec, 'mode', 4, 'dE', fwhmESplit, 'qlabel', {'M' '\Gamma' 'M'})
ylim([1 7]);
clim([-1 1])
legend off;
title('')
ax = gca;
ax.FontSize = 18;
axis square
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
set(gca, 'FontSize', 18)
colormap(redblue(255))

spec = MnF2Splitting.spinwave({[-0.5,0.5,1.5] [0,0,1.5] [0.5,0.5,1.5] 223}, 'formfact', true, 'hermit', true);
spec = sw_neutron(spec, 'pol', true, 'uv', {[-0.5 0.5 1.5] [0.5 0.5 1.5]});
spec = sw_egrid(spec, 'Evect', linspace(0, 10, 1000), 'component', 'Pxy+Pxz', 'T', 0, 'ImagChk', true);

nexttile()
sw_plotspec(spec, 'mode', 4, 'dE', fwhmESplit/3, 'qlabel', {'A' 'Z' 'A'})
ylim([6 8]);
clim([-1 1])
legend off;
title('')
ax = gca;
ax.FontSize = 18;
axis square
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
set(gca, 'FontSize', 18)
colormap(redblue(255))

set(f, 'Position', [2, 3.3, 15, 13.7], 'Units', 'centimeters')

disp(['Max delta E: ', num2str(max(spec.omega(1,:) - spec.omega(2,:)))])

saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig1/Splitting.pdf'))

%% Figure 1: Cut to determine when we get 0.25 meV splitting

J1 = fitStr.x(1); J2 = fitStr.x(2); J3 = fitStr.x(3); J7 = 0.006; J8 = -0.006; Dc = fitStr.x(4); fwhmESplit = 0.2;
%J1 = fitStr.x(1); J2 = fitStr.x(2); J3 = fitStr.x(3); J7 = 0.005; J8 = -0.005; J9 = -0.005; J10 = 0.005; Dc = fitStr.x(4); fwhmESplit = 0.4;

MnF2Splitting = spinw;
MnF2Splitting.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled', [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Splitting.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
MnF2Splitting.gencoupling('maxDistance', 10);
MnF2Splitting.addmatrix('label', 'J1', 'value', J1);
MnF2Splitting.addmatrix('label', 'J2', 'value', J2);
MnF2Splitting.addmatrix('label', 'J3', 'value', J3);
MnF2Splitting.addmatrix('label', 'J7', 'value', J7);
MnF2Splitting.addmatrix('label', 'J8', 'value', J8);
% MnF2Splitting.addmatrix('label', 'J9', 'value', J9);
% MnF2Splitting.addmatrix('label', 'J10', 'value', J10);
MnF2Splitting.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
MnF2Splitting.addcoupling('mat', 'J1', 'bond', 1);
MnF2Splitting.addcoupling('mat', 'J2', 'bond', 2);
MnF2Splitting.addcoupling('mat', 'J3', 'bond', 3);
MnF2Splitting.addcoupling('mat', 'J7', 'bond', 7);
MnF2Splitting.addcoupling('mat', 'J8', 'bond', 8);
% MnF2Splitting.addcoupling('mat', 'J9', 'bond', 9);
% MnF2Splitting.addcoupling('mat', 'J10', 'bond', 10);
MnF2Splitting.addaniso('Dc', 1);
MnF2Splitting.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -5/2 5/2]);

% MnF2Splitting.table('bond',1:10)

spec = MnF2Splitting.spinwave({[1,1,0.5] [0.5,0.5,0.5] 223}, 'formfact', true, 'hermit', true); % Update indSplit if change binning.
spec = sw_neutron(spec);
spec = sw_egrid(spec, 'Evect', linspace(0, 10, 5000), 'component', 'SPerp', 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmESplit);

figure
sw_plotspec(spec, 'mode', 4, 'qlabel', {'A' 'Z'})
ylim([6.5 7.5]);
clim([0 2])
legend off;
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
colormap('turbo')

indSplit = 112; % Index for 0.75,0.75,0.5. Update if change binning.
ESim = ((spec.Evect(2:end)-spec.Evect(1:end-1))/2 + spec.Evect(1:end-1))'; % Convert from bin edge to bin center
intSim = spec.swConv(:,indSplit);

f = fit(ESim, intSim, 'gauss1', 'StartPoint', [0.08, 6.9, 0.15]);

figure
hold on
plot(f, ESim, intSim)
ylabel('\itI\rm (arb. units)')
xlabel('\rm\Delta\itE\rm (meV)')
xlim([6.5 7.5]);
box on
hold off

disp(['FWHM: ', num2str(f.c1*2*sqrt(2*log(2))/sqrt(2))]) % MATLAB's gauss1 doesn't have the factor of 1/2 in exponential so need to divide sigma by sqrt(2). Can confirm by setting fwhm=1 in convolution calculation.
disp(['Max delta E: ', num2str(max(spec.omega(1,:) - spec.omega(2,:)))])

%% Figure 1: Calculate splitting at (3/4,3/4,1/2) analytically from https://arxiv.org/pdf/2410.10771v1. See 01.11.2024 notes.

S = 5/2;
J7a = 0.006;
J7b = -0.006;
a = 4.8736;
kx = 2*pi/a*3/4; ky = 2*pi/a*3/4;
deltaE = abs(4*S*(J7a-J7b)*sin(a*kx)*sin(a*ky));

disp(['Splitting: ' num2str(deltaE) ' meV'])

%% Figure 2: Plot the simulated high-symmetry paths (tutorial 35)

spec = MnF2Fit.spinwave({[0,0,1] [0.5,0.5,1] [0.5,0.5,0.5] [1,1,0.5] [1,1,0] 223}, 'formfact', true,... 
    'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(0, 9, 1000), 'component',...
    'SPerp', 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmE);
figure
sw_plotspec(spec, 'mode', 4, 'qlabel', {'\Gamma' 'M' 'A' 'Z' '\Gamma'})
ylim([0 8.4]);
clim([0 5.5]);
colormap(gca, 'turbo');
legend off;
hold on
% Gamma to M data
wind1 = [5 57:65]; % Indices from datafile
Qabs1 = sqrt(sum((bsxfun(@minus,Q(wind1,:),[0 0 1])*MnF2Fit.rl).^2,2)); % Distance starting from 0
dist1 = sqrt(sum(([0.5 0.5 0]*MnF2Fit.rl).^2)); % Total distance by end
% M to A data
wind2 = flip(16:26);
Qabs2 = sqrt(sum((bsxfun(@minus,Q(wind2,:),[0.5 0.5 1])*MnF2Fit.rl).^2,2)) + dist1;
dist2 = dist1 + sqrt(sum(([0 0 0.5]*MnF2Fit.rl).^2));
% A to Z data
wind3 = 48:56;
Qabs3 = sqrt(sum((bsxfun(@minus,Q(wind3,:),[0.5 0.5 0.5])*MnF2Fit.rl).^2,2)) + dist2;
dist3 = dist2 + sqrt(sum(([0.5 0.5 0]*MnF2Fit.rl).^2));
% Z to Gamma data
wind4 = flip(27:35);
Qabs4 = sqrt(sum((bsxfun(@minus,Q(wind4,:),[1 1 0.5])*MnF2Fit.rl).^2,2)) + dist3;
wabs = w([wind1, wind2, wind3, wind4])';
wabserr = werr([wind1, wind2, wind3, wind4])';
Qabs = [Qabs1; Qabs2; Qabs3; Qabs4];
intAbs = inten([wind1, wind2, wind3, wind4])';
for i=1:length(wabs) % Loop through to vary marker size based on intensity
    plot3(Qabs(i), wabs(i), 1e2+wabs(i)*0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 2/max(intAbs)*intAbs(i)+4)
end
plot3([Qabs, Qabs]', [-wabserr, wabserr]'+wabs', [1e2+wabs*0, 1e2+wabs*0]', 'w-', 'MarkerFaceColor', 'w', 'LineWidth', 2)
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 10, 5]);
for i = 1:length(ax.XTick) % Add vertical lines at each Q that are solid white lines
    xline(ax.XTick(i), 'Color', 'w', 'LineWidth', 1.5)
end
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/HighSymHHL.pdf'))

spec = MnF2Fit.spinwave({[1,0,0] [1.5,0,0] [1.5,0,0.5] [1,0,0.5] 223}, 'formfact', true,... 
    'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(0, 9, 1000), 'component',...
    'SPerp', 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmE);
figure
sw_plotspec(spec, 'mode', 4, 'qlabel', {'\Gamma' 'X' 'R' 'Z'})
ylim([0 8.4]);
clim([0 4]);
colormap(gca, 'turbo');
legend off;
hold on
% Gamma to X data
wind1 = 74:84; % Indices from datafile
Qabs1 = sqrt(sum((bsxfun(@minus,Q(wind1,:),[1 0 0])*MnF2Fit.rl).^2,2)); % Distance starting from 0
dist1 = sqrt(sum(([0.5 0 0]*MnF2Fit.rl).^2)); % Total distance by end
% X to R data
wind2 = 113:122;
Qabs2 = sqrt(sum((bsxfun(@minus,Q(wind2,:),[1.5 0 0])*MnF2Fit.rl).^2,2)) + dist1;
dist2 = dist1 + sqrt(sum(([0 0 0.5]*MnF2Fit.rl).^2));
% R to Z data
wind3 = flip(91:101);
Qabs3 = sqrt(sum((bsxfun(@minus,Q(wind3,:),[1.5 0 0.5])*MnF2Fit.rl).^2,2)) + dist2;
wabs = w([wind1, wind2, wind3])';
wabserr = werr([wind1, wind2, wind3])';
Qabs = [Qabs1; Qabs2; Qabs3];
intAbs = inten([wind1, wind2, wind3])';
for i=1:length(wabs) % Loop through to vary marker size based on intensity
    plot3(Qabs(i), wabs(i), 1e2+wabs(i)*0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 2/max(intAbs)*intAbs(i)+4)
end
plot3([Qabs, Qabs]', [-wabserr, wabserr]'+wabs', [1e2+wabs*0, 1e2+wabs*0]', 'w-', 'MarkerFaceColor', 'w', 'LineWidth', 2)
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 7.5, 5]);
set(gca, 'YTickLabel', [])
for i = 1:length(ax.XTick) % Add vertical lines at each Q that are solid white lines
    xline(ax.XTick(i), 'Color', 'w', 'LineWidth', 1.5)
end
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/HighSymH0L.pdf'))

spec = MnF2Fit.spinwave({[1,0.5,0] [1.5,0.5,0] 223}, 'formfact', true,... 
    'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(0, 9, 1000), 'component',...
    'SPerp', 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmE);
figure
sw_plotspec(spec, 'mode', 4, 'qlabel', {'X' 'M'})
ylim([0 8.75]);
clim([0 1]);
colormap(gca, 'turbo');
legend off;
hold on
% X to M data. One path so plots against H rather than Qabs.
wind = 132:142;
intAbs = inten(wind)';
for i=1:length(wind) % Loop through to vary marker size based on intensity
    plot3(Q(wind(i),1)-1, w(wind(i)), 1e2+w(wind(i))*0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 2/max(intAbs)*intAbs(i)+4)
end
plot3([Q(wind,1)-1, Q(wind,1)-1]', [-werr(wind); werr(wind)]+w(wind), [1e2+w(wind)*0; 1e2+w(wind)*0], 'w-', 'MarkerFaceColor', 'w', 'LineWidth', 2)
xticks([0 0.5])
xticklabels({'X' 'M'}) % qlabel fails for two points
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 2.6, 5]);
set(gca, 'YTickLabel', [])
for i = 1:length(ax.XTick) % Add vertical lines at each Q that are solid white lines
    xline(ax.XTick(i), 'Color', 'w', 'LineWidth', 1.5)
end
disp(['Max delta E: ', num2str(max(spec.omega(1,:) - spec.omega(2,:)))])
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/HighSymHK0.pdf'))

%% Figure 2: Plot the simulated constant-E slices (tutorial 10 without Horace, 21.10.2024 notes)

maskLim1 = [0.5 1.60]; % How far to mask simulation from Gamma
maskLim2 = [0.55 2.12]; % How far to mask simulation from Gamma
maskLim3 = [0.7 2.05]; % How far to mask simulation from Gamma
Ecut1 = [3.95 4.05]; % meV
Ecut2 = [4.95 5.05]; % meV
Ecut3 = [5.95 6.05]; % meV
dir1 = 'U:/Data/Figures/MnF2/Manuscript/Fig2/HK04meV.csv';
dir2 = 'U:/Data/Figures/MnF2/Manuscript/Fig2/HK05meV.csv';
dir3 = 'U:/Data/Figures/MnF2/Manuscript/Fig2/HK06meV.csv';
scale = 8e-2; % Overall scale factor for simulation
cMax = 1; % Scale factor for data and simulation
modturbo = [1 1 1; turbo]; % Make negative/missing values white.

nQ = 301;
nE = 1001;
Qhv = linspace(0,1.7,nQ);
Qkv = linspace(0,1.7,nQ);
Qlv = 0;
[Qk, Qh, Ql] = ndgrid(Qkv, Qhv, Qlv);
Q = [Qh(:) Qk(:) Ql(:)]';

spec = MnF2Fit.spinwave(Q, 'formfact', true, 'hermit', true);
Ev = linspace(3,8,nE);
spec = sw_egrid(spec, 'component', 'SPerp', 'Evect', Ev, 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmE);

spec3D = reshape(spec.swConv,nE-1,nQ,nQ);

Eidx1 = find(Ev>Ecut1(1) & Ev<Ecut1(2));
int1 = squeeze(sum(spec3D(Eidx1,:,:),1))/numel(Eidx1)/(Ev(2)-Ev(1));
Eidx2 = find(Ev>Ecut2(1) & Ev<Ecut2(2));
int2 = squeeze(sum(spec3D(Eidx2,:,:),1))/numel(Eidx2)/(Ev(2)-Ev(1));
Eidx3 = find(Ev>Ecut3(1) & Ev<Ecut3(2));
int3 = squeeze(sum(spec3D(Eidx3,:,:),1))/numel(Eidx3)/(Ev(2)-Ev(1));

% Plot the simulated constant-E slices with the data from Python. Data is
% already binned for plotting.
file1 = readmatrix(dir1, 'NumHeaderLines', 1);
QhDatRaw1 = file1(:,4);
QkDatRaw1 = file1(:,5);
intDatRaw1 = file1(:,12);
QhUn1 = length(uniquetol(QkDatRaw1)); QkUn1 = length(uniquetol(QhDatRaw1));
QhDat1 = reshape(QhDatRaw1, QhUn1, QkUn1);
QkDat1 = reshape(QkDatRaw1, QhUn1, QkUn1);
intDat1 = reshape(intDatRaw1, QhUn1, QkUn1);
% Crop data so square matrix
QhDat1 = QhDat1(1:min([QhUn1, QkUn1]), 1:min([QhUn1, QkUn1]));
QkDat1 = QkDat1(1:min([QhUn1, QkUn1]), 1:min([QhUn1, QkUn1]));
intDat1 = intDat1(1:min([QhUn1, QkUn1]), 1:min([QhUn1, QkUn1]));
% Plot the data and the simulation
figure
hold on
intPlt1 = tril(int1*scale, -1) - triu(ones(size(int1))); % Mask pixels in other region
intPlt1(sqrt(Qh.^2+Qk.^2)<maskLim1(1) | sqrt(Qh.^2+Qk.^2)>maskLim1(2))=-1; % Mask pixels in region beyond data
p1 = surf(Qh, Qk, intPlt1);
p1.EdgeColor = 'none';
QhDatPlt1 = QhDat1; QkDatPlt1 = QkDat1; intDatPlt1 = triu(intDat1, 1) - tril(ones(size(intDat1))); % Mask pixels in other region
p2 = surf(QhDatPlt1, QkDatPlt1, intDatPlt1);
p2.EdgeColor = 'none';
xlabel('(\it{H}\rm{00) (rlu)}')
ylabel('(0\it{K}\rm{0) (rlu)}')
xlim([0 1.5])
ylim([0 1.5])
clim([-1e-2 cMax]) % Place the values super close to 0 e.g. from the calculation just above the first pixel of the colormap which is white
colormap(gca, modturbo)
c = colorbar;
view(2)
ax = gca;
ax.FontSize = 16;
c.Ticks = [0 cMax];
c.TickLabels = {0 1};
axis square
box on
set(gca, 'layer', 'top')
plot3([0.5, 0.5], [0, 1.5], [1e3, 1e3], '--w', 'LineWidth', 2)
plot3([1, 1], [0, 1.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 1.5], [0.5, 0.5], [1e3, 1e3], '--w', 'LineWidth', 2)
plot3([0, 1.5], [1, 1], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 1], [1, 0], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0.5, 1.5], [1.5, 0.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 0.5], [1, 1.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([1, 1.5], [0, 0.5], [1e3, 1e3], '--w', 'LineWidth', 1)
text(0.01, 1.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
text(0.51, 1.07, 'X', 'Color', 'w', 'FontSize', 16);
text(0.31, 0.56, 'M', 'Color', 'w', 'FontSize', 16);
text(1.01, 0.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
text(1.01, 0.57, 'X', 'Color', 'w', 'FontSize', 16);
text(1.01, 1.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/ESliceHK0', num2str(round(mean(Ecut1), 1)), 'meV.pdf'))
hold off

% Plot the simulated constant-E slices with the data from Python. Data is
% already binned for plotting.
file2 = readmatrix(dir2, 'NumHeaderLines', 1);
QhDatRaw2 = file2(:,4);
QkDatRaw2 = file2(:,5);
intDatRaw2 = file2(:,12);
QhUn2 = length(uniquetol(QkDatRaw2)); QkUn2 = length(uniquetol(QhDatRaw2));
QhDat2 = reshape(QhDatRaw2, QhUn2, QkUn2);
QkDat2 = reshape(QkDatRaw2, QhUn2, QkUn2);
intDat2 = reshape(intDatRaw2, QhUn2, QkUn2);
% Crop data so square matrix
QhDat2 = QhDat2(1:min([QhUn2, QkUn2]), 1:min([QhUn2, QkUn2]));
QkDat2 = QkDat2(1:min([QhUn2, QkUn2]), 1:min([QhUn2, QkUn2]));
intDat2 = intDat2(1:min([QhUn2, QkUn2]), 1:min([QhUn2, QkUn2]));
% Plot the data and the simulation
figure
hold on
intPlt2 = tril(int2*scale, -1) - triu(ones(size(int2))); % Mask pixels in other region
intPlt2(sqrt(Qh.^2+Qk.^2)<maskLim2(1) | sqrt(Qh.^2+Qk.^2)>maskLim2(2))=-1; % Mask pixels in region beyond data
p1 = surf(Qh, Qk, intPlt2);
p1.EdgeColor = 'none';
QhDatPlt2 = QhDat2; QkDatPlt2 = QkDat2; intDatPlt2 = triu(intDat2, 1) - tril(ones(size(intDat2))); % Mask pixels in other region
p2 = surf(QhDatPlt2, QkDatPlt2, intDatPlt2);
p2.EdgeColor = 'none';
xlabel('(\it{H}\rm{00) (rlu)}')
ylabel('(0\it{K}\rm{0) (rlu)}')
xlim([0 1.5])
ylim([0 1.5])
clim([-1e-2 cMax])
colormap(gca, modturbo)
c = colorbar;
view(2)
ax = gca;
ax.FontSize = 16;
c.Ticks = [0 cMax];
c.TickLabels = {0 1};
axis square
box on
set(gca, 'layer', 'top')
plot3([0.5, 0.5], [0, 1.5], [1e3, 1e3], '--w', 'LineWidth', 2)
plot3([1, 1], [0, 1.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 1.5], [0.5, 0.5], [1e3, 1e3], '--w', 'LineWidth', 2)
plot3([0, 1.5], [1, 1], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 1], [1, 0], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0.5, 1.5], [1.5, 0.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 0.5], [1, 1.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([1, 1.5], [0, 0.5], [1e3, 1e3], '--w', 'LineWidth', 1)
text(0.01, 1.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
text(0.51, 1.07, 'X', 'Color', 'w', 'FontSize', 16);
text(0.31, 0.56, 'M', 'Color', 'w', 'FontSize', 16);
text(1.01, 0.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
text(1.01, 0.57, 'X', 'Color', 'w', 'FontSize', 16);
text(1.01, 1.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/ESliceHK0', num2str(round(mean(Ecut2), 1)), 'meV.pdf'))
hold off

% Plot the simulated constant-E slices with the data from Python. Data is
% already binned for plotting.
file3 = readmatrix(dir3, 'NumHeaderLines', 1);
QhDatRaw3 = file3(:,4);
QkDatRaw3 = file3(:,5);
intDatRaw3 = file3(:,12);
QhUn3 = length(uniquetol(QkDatRaw3)); QkUn3 = length(uniquetol(QhDatRaw3));
QhDat3 = reshape(QhDatRaw3, QhUn3, QkUn3);
QkDat3 = reshape(QkDatRaw3, QhUn3, QkUn3);
intDat3 = reshape(intDatRaw3, QhUn3, QkUn3);
% Crop data so square matrix
QhDat3 = QhDat3(1:min([QhUn3, QkUn3]), 1:min([QhUn3, QkUn3]));
QkDat3 = QkDat3(1:min([QhUn3, QkUn3]), 1:min([QhUn3, QkUn3]));
intDat3 = intDat3(1:min([QhUn3, QkUn3]), 1:min([QhUn3, QkUn3]));
% Plot the data and the simulation
figure
hold on
intPlt3 = tril(int3*scale, -1) - triu(ones(size(int3))); % Mask pixels in other region
intPlt3(sqrt(Qh.^2+Qk.^2)<maskLim3(1) | sqrt(Qh.^2+Qk.^2)>maskLim3(2))=-1; % Mask pixels in region beyond data
p1 = surf(Qh, Qk, intPlt3);
p1.EdgeColor = 'none';
QhDatPlt3 = QhDat3; QkDatPlt3 = QkDat3; intDatPlt3 = triu(intDat3, 1) - tril(ones(size(intDat3))); % Mask pixels in other region
p2 = surf(QhDatPlt3, QkDatPlt3, intDatPlt3);
p2.EdgeColor = 'none';
grid off
xlabel('(\it{H}\rm{00) (rlu)}')
ylabel('(0\it{K}\rm{0) (rlu)}')
xlim([0 1.5])
ylim([0 1.5])
clim([-1e-2 cMax])
colormap(gca, modturbo)
c = colorbar;
view(2)
ax = gca;
ax.FontSize = 16;
c.Ticks = [0 cMax];
c.TickLabels = {0 1};
ylabel(c, '\itI\rm (arb. units)')
axis square
box on
set(gca, 'layer', 'top')
plot3([0.5, 0.5], [0, 1.5], [1e3, 1e3], '--w', 'LineWidth', 2)
plot3([1, 1], [0, 1.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 1.5], [0.5, 0.5], [1e3, 1e3], '--w', 'LineWidth', 2)
plot3([0, 1.5], [1, 1], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 1], [1, 0], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0.5, 1.5], [1.5, 0.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([0, 0.5], [1, 1.5], [1e3, 1e3], '--w', 'LineWidth', 1)
plot3([1, 1.5], [0, 0.5], [1e3, 1e3], '--w', 'LineWidth', 1)
text(0.01, 1.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
text(0.51, 1.07, 'X', 'Color', 'w', 'FontSize', 16);
text(0.31, 0.56, 'M', 'Color', 'w', 'FontSize', 16);
text(1.01, 0.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
text(1.01, 0.57, 'X', 'Color', 'w', 'FontSize', 16);
text(1.01, 1.1, '\Gamma', 'Color', 'w', 'FontSize', 16);
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/ESliceHK0', num2str(round(mean(Ecut3), 1)), 'meV.pdf'))
hold off

%% Figure 3: Export dispersion for high-resolution 0 T slice

% Fitted model
spec = MnF2Fit.spinwave({[0.45,0.45,0.5] [1.05,1.05,0.5] 223}, 'formfact', true, 'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(6.0, 7.5, 1000), 'component', 'SPerp', 'T', 1.6, 'ImagChk', true);
omegaFit = spec.omega(1:2,:);
save('U:\Data\MATLAB\MnF2\omegaFit.txt', 'omegaFit', '-ascii')
hklFit = spec.hkl;
save('U:\Data\MATLAB\MnF2\hklFit.txt', 'hklFit', '-ascii')

%% Figure 3: Split model
J1 = fitStr.x(1); J2 = fitStr.x(2); J3 = fitStr.x(3); J4 = 0;...
    J7 = 0.02; J8 = -0.02; J9 = -0.00; J10 = 0.00; J7B = -0.00;...
    J8B = 0.00; J9B = -0.00; J10B = 0.00; Dc = fitStr.x(4);...
    fwhmSplit = 0.2;

MnF2Split = spinw;
MnF2Split.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled', [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Split.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});

MnF2Split.gencoupling('maxDistance', 10);
MnF2Split.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
MnF2Split.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
MnF2Split.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
MnF2Split.addmatrix('label', 'J7', 'value', J7, 'color', 'b');
MnF2Split.addmatrix('label', 'J8', 'value', J8, 'color', 'g');
MnF2Split.addmatrix('label', 'J9', 'value', J9);
MnF2Split.addmatrix('label', 'J10', 'value', J10);
MnF2Split.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))
MnF2Split.addcoupling('mat', 'J1', 'bond', 1);
MnF2Split.addcoupling('mat', 'J2', 'bond', 2);
MnF2Split.addcoupling('mat', 'J3', 'bond', 3);
MnF2Split.addcoupling('mat', 'J7', 'bond', 7);
MnF2Split.addcoupling('mat', 'J8', 'bond', 8);
MnF2Split.addcoupling('mat', 'J9', 'bond', 9);
MnF2Split.addcoupling('mat', 'J10', 'bond', 10);
MnF2Split.addaniso('Dc', 1);
MnF2Split.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);

spec = MnF2Split.spinwave({[0.45,0.45,0.5] [1.05,1.05,0.5] 223}, 'formfact', true, 'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(6.0, 7.5, 1000), 'component', 'SPerp', 'T', 1.6, 'ImagChk', true);
omegaSplit = spec.omega(1:2,:);
save('U:\Data\MATLAB\MnF2\omegaSplit.txt', 'omegaSplit', '-ascii')
hklSplit = spec.hkl;
save('U:\Data\MATLAB\MnF2\hklSplit.txt', 'hklSplit', '-ascii')

% figure
% subplot(2, 1, 1)
% sw_plotspec(spec, 'mode', 1, 'axLim', [0 10], 'qlabel', {'\Gamma' 'X' 'M' '\Gamma' 'Z' 'R' 'A' 'Z'})
% %sw_plotspec(spec, 'mode', 1, 'axLim', [0 7])
% subplot(2, 1, 2)
% sw_plotspec(spec, 'mode', 2, 'axLim', [0 10], 'qlabel', {'\Gamma' 'X' 'M' '\Gamma' 'Z' 'R' 'A' 'Z'})
% %sw_plotspec(spec, 'mode', 2, 'axLim', [0 7])

figure
sw_plotspec(spec, 'mode', 4, 'dE', fwhmSplit)
ylim([6.4 7.2]);
clim([0 3])
colormap(gca, 'turbo');
legend off;
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')

%% Figure 4: Plot the simulated 10 T spectrum

% Adjust the g-tensor
% MnF2Fit.addmatrix('label', 'g', 'value', diag([2.2,2.2,2.2]));
% MnF2Fit.gencoupling;
% MnF2Fit.addg('g', 'MMn2');

% Apply a magnetic field along (H,-H,L)
MnF2Fit.field([10,-10,0]/sqrt(2));

% Relax magnetic structure
MnF2Fit.optmagsteep('nRun', 1e4)

spec10T = MnF2Fit.spinwave({[0,0,1] [0.5,0.5,1] [0.5,0.5,0.5] [1,1,0.5]...
    [1,1,0] 223}, 'formfact', true, 'hermit', true, 'gtensor', false);
spec10T = sw_egrid(spec10T, 'Evect', linspace(0, 9, 1000), 'component',...
    'SPerp', 'T', 1.6, 'ImagChk', true);
spec10T = sw_instrument(spec10T, 'dE', fwhmE);
figure
sw_plotspec(spec10T, 'mode', 4, 'qlabel', {'\Gamma' 'M' 'A' 'Z' '\Gamma'})
ylim([0 8.4]);
clim([0 5]);
colormap(gca, 'turbo');
c = colorbar;
% c.Ticks = [0 5];
% c.TickLabels = {0 1};
ylabel(c, '\itI\rm (arb. units)')
legend off;
hold on
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 10, 5]);
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig4/HighSym10T.pdf'))

% Figure 4: Export dispersion for simulated 10 T spectrum
omegaFit = spec10T.omega(1:2,:);
save('U:\Data\MATLAB\MnF2\omegaFit10T.txt', 'omegaFit', '-ascii')
hklFit = spec10T.hkl;
save('U:\Data\MATLAB\MnF2\hklFit10T.txt', 'hklFit', '-ascii')

%% Figure 4: Cut through the simulated 10 T spectrum, export for Python

% Include dQ FWHM
%fwhmQ = 0.0622; % inverse Angstroms determined from MJOLNIR calculation
fwhmQ = 0.0; % inverse Angstroms

spec10T = MnF2Fit.spinwave({[0,0,1] [0.5,0.5,1] 223}, 'formfact', true, 'hermit', true, 'gtensor', false);
%spec10T = MnF2Fit.spinwave({[0,0,1] [0,0,1.5] 223}, 'formfact', true, 'hermit', true, 'gtensor', false); % Walk radially to use SpinW Q-convolution
spec10T = sw_egrid(spec10T, 'Evect', linspace(0, 9, 1000), 'component', 'SPerp', 'T', 1.6, 'ImagChk', true);
spec10T = sw_instrument(spec10T, 'dE', fwhmE);
%spec10T = sw_instrument(spec10T, 'dE', fwhmE, 'dQ', fwhmQ);

figure
sw_plotspec(spec10T, 'mode', 4)
ylim([0 8.4]);
clim([0 5]);
colormap(gca, 'turbo');
legend off;
hold on
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 10, 5]);

ESim = (spec10T.Evect(2:end)-spec10T.Evect(1:end-1))/2 + spec10T.Evect(1:end-1); % Convert from bin edge to bin center
intSim = spec10T.swConv(:,1); % Intensity at Gamma
%intSim = sum(spec10T.swConv(:,1:10), 2); % Summed intensity going out to half the total bin width from the E cuts in Python

figure
hold on
plot(ESim, intSim)
ylabel('\itI\rm (arb. units)')
xlabel('\rm\Delta\itE\rm (meV)')
xlim([0.5, 2.5])
box on
hold off

save('U:\Data\MATLAB\MnF2\ESim10T.txt', 'ESim', '-ascii')
save('U:\Data\MATLAB\MnF2\intSim10T.txt', 'intSim', '-ascii')

%% Fit J1, J2, J3 with dipole-dipole interaction rather than single-ion

% Global Fit: Develop SpinW model for Hamiltonian

clear, close all

intDist = 20; % Max coupling distance in Angstroms, used for dipole-dipole. Iterate for convergence.
% J1 = 2*-0.028; J2 = 2*0.152; J3 = 2*0.004; % Nikotin 1969
J1 = 2*-0.0354; J2 = 2*0.1499; J3 = 2*0.000; % CAMEA 2023
% J1 = -0.064444; J2 = 0.31081; J3 = -0.0020532; % Previous values
% J1 = -0.08; J2 = 0.32; J3 = -0.006; % Previous values

% Instrumental E-resolution extracted from MJOLNIR for HHL along (001).
% Will be fit to polynomial in sw_instrument.
sigE = [0.3090, 0.0947; 1.0043, mean([0.0917, 0.0948]);...
    2.0002, mean([0.0889, 0.0919]); 2.9962, mean([0.1340, 0.1372]);...
    3.9921, mean([0.1321, 0.1351]); 5.0069, mean([0.1793, 0.1827]);...
    6.0028, mean([0.1779, 0.1811]); 6.9988, mean([0.2276, 0.2310]);...
    7.9947, mean([0.2264, 0.2298])];
fwhmE = sigE; fwhmE(:,2) = sigE(:,2)*2*sqrt(2*log(2));

MnF2Fit = spinw;
MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
    [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label',...
    {'MMn2' 'F'}, 'color', {'red' 'yellow'});

MnF2Fit.gencoupling('maxDistance', intDist); % Larger distance for dipole-dipole convergence
% MnF2Fit.table('bond',1:10)

MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');

MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);

MnF2Fit.coupling.rdip = intDist; % Dipole-dipole interaction
MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);

% Relax magnetic structure
% MnF2Fit.optmagsteep('nRun', 1e4)
% MnF2Fit.table('mag')

% plot(MnF2Fit, 'range', [1 1 1], 'bondMode', 'line', 'bondLinewidth0', 3)

%% Global Fit: Prepare fitting to the dispersion (tutorial 35)

x1 = [J1 J2 J3];
% x1 = [J1 J2];
MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2' 'J3'}, 'init', true)
% MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2'}, 'init', true)

MnF2Fit.table('matrix')

pref = swpref;
pref.tid = 0;

par_fit          = struct;
par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';

par_fit.Evect    = linspace(0, 10, 1001);
% par_fit.Evect    = linspace(0, 10, 81);
par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                    'mat', {'J1' 'J2' 'J3'}, 'init', true);
par_fit.xmin      = [-0.6 -0.6 -0.2];
par_fit.xmax      = [0.6 0.6 0.2];
% par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
%                     'mat', {'J1' 'J2'}, 'init', true);
% par_fit.xmin      = [-0.6 -0.6];
% par_fit.xmax      = [0.6 0.6];
par_fit.x0        = x1;
par_fit.plot      = true;
par_fit.hermit    = false;
par_fit.optimizer = 'pso';
par_fit.maxiter   = 1e2;
par_fit.nrun      = 1;
par_fit.TolX      = 10^(-4);
par_fit.TolFun    = 10^(-5);

%% Global Fit: Fit to the dispersion (tutorial 35)

fitStr = MnF2Fit.fitspec(par_fit);

disp(['p: ' num2str(fitStr.x)])
disp(['chi2: ' num2str(fitStr.redX2)])

%% Global Fit: Extract uncertainties

% J1
delJ1 = linspace(-0.01, 0.01, 20);
J1chi2 = zeros(size(delJ1));
J1Iter = zeros(size(delJ1));
T = sw_readtable('U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt');
dof = length(T) - 2 + 1; % https://spinw.org/spinwdoc/spinw_fitspec
for i=1:length(delJ1)
    J1 = fitStr.x(1) + delJ1(i); J2 = fitStr.x(2);...
        J3 = fitStr.x(3);
    
    MnF2Fit = spinw;
    MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
        [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
    MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0],...
        'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
    MnF2Fit.gencoupling('maxDistance', intDist);
    MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
    MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
    MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
    MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
    MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
    MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
    MnF2Fit.coupling.rdip = intDist; % Dipole-dipole interaction
    MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);
    
    x1 = [J2 J3];
    MnF2Fit.matparser('param', x1, 'mat', {'J2' 'J3'}, 'init', true)
    
    MnF2Fit.table('matrix')
    
    pref = swpref;
    pref.tid = 0;
    
    par_fit          = struct;
    par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';
    
    par_fit.Evect    = linspace(0, 10, 501);
    par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                        'mat', {'J2' 'J3'}, 'init', true);
    par_fit.xmin      = [-0.6 -0.2];
    par_fit.xmax      = [0.6 0.2];
    par_fit.x0        = x1;
    par_fit.plot      = false;
    par_fit.hermit    = false;
    par_fit.optimizer = 'pso';
    par_fit.maxiter   = 1e2;
    par_fit.nrun      = 1;
    par_fit.TolX      = 10^(-4);
    par_fit.TolFun    = 10^(-5);
    
    fitStrIter = MnF2Fit.fitspec(par_fit);

    J1Iter(i) = fitStr.x(1) + delJ1(i);
    J1chi2(i) = fitStrIter.redX2;
    
    disp(['i: ' num2str(i)])
    disp(['p: ' num2str(fitStrIter.x)])
    disp(['chi2R: ' num2str(fitStrIter.redX2)])
end
chi2RThr = (1 + 1/dof)*min(J1chi2); % Threshold reduced chi-squared
indClose = find(J1chi2<chi2RThr); % J1's within the threshold
sigJ1 = (J1Iter(indClose(end)+1) - J1Iter(indClose(1)-1))/2; % Go to points just beyond intercept with threshold, subtract, divide by two to slightly overestimate uncertainty
figure
hold on
scatter(J1Iter, J1chi2)
yline(chi2RThr)
xlabel('J1 (meV)')
ylabel('Chi2R')
title(['J1: ' num2str(fitStr.x(1)) '. sigJ1: ' num2str(sigJ1)])
box on
hold off
saveas(gcf, strcat('U:/Data/Figures/MnF2/Fitting/UncJ1Dip.pdf'))
disp(['J1: ' num2str(fitStr.x(1))])
disp(['sigJ1: ' num2str(sigJ1)])

% J2
delJ2 = linspace(-0.01, 0.01, 20);
J2chi2 = zeros(size(delJ2));
J2Iter = zeros(size(delJ2));
for i=1:length(delJ2)
    J1 = fitStr.x(1); J2 = fitStr.x(2) + delJ2(i);...
        J3 = fitStr.x(3);
    
    MnF2Fit = spinw;
    MnF2Fit.genlattice('lat_const', [4.87380 4.87380 3.31070], 'angled',...
        [90 90 90], 'sym', 136);
    MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0],...
        'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
    MnF2Fit.gencoupling('maxDistance', intDist);
    MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
    MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
    MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
    MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
    MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
    MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
    MnF2Fit.coupling.rdip = intDist; % Dipole-dipole interaction
    MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);
    
    x1 = [J1 J3];
    MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J3'}, 'init', true)
    
    MnF2Fit.table('matrix')
    
    pref = swpref;
    pref.tid = 0;
    
    par_fit          = struct;
    par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';
    
    par_fit.Evect    = linspace(0, 10, 501);
    par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                        'mat', {'J1' 'J3' 'Dc(3,3)'}, 'init', true);
    par_fit.xmin      = [-0.6 -0.2];
    par_fit.xmax      = [0.6 0.2];
    par_fit.x0        = x1;
    par_fit.plot      = false;
    par_fit.hermit    = false;
    par_fit.optimizer = 'pso';
    par_fit.maxiter   = 1e2;
    par_fit.nrun      = 1;
    par_fit.TolX      = 10^(-4);
    par_fit.TolFun    = 10^(-5);
    
    fitStrIter = MnF2Fit.fitspec(par_fit);

    J2Iter(i) = fitStr.x(2) + delJ2(i);
    J2chi2(i) = fitStrIter.redX2;
    
    disp(['i: ' num2str(i)])
    disp(['p: ' num2str(fitStrIter.x)])
    disp(['chi2R: ' num2str(fitStrIter.redX2)])
end
chi2RThr = (1 + 1/dof)*min(J2chi2); % Threshold reduced chi-squared
indClose = find(J2chi2<chi2RThr); % J2's within the threshold
sigJ2 = (J2Iter(indClose(end)+1) - J2Iter(indClose(1)-1))/2; % Go to points just beyond intercept with threshold, subtract, divide by two to slightly overestimate uncertainty
figure
hold on
scatter(J2Iter, J2chi2)
yline(chi2RThr)
xlabel('J2 (meV)')
ylabel('Chi2R')
title(['J2: ' num2str(fitStr.x(2)) '. sigJ2: ' num2str(sigJ2)])
box on
hold off
saveas(gcf, strcat('U:/Data/Figures/MnF2/Fitting/UncJ2Dip.pdf'))
disp(['J2: ' num2str(fitStr.x(2))])
disp(['sigJ2: ' num2str(sigJ2)])

% J3
delJ3 = linspace(-0.01, 0.01, 20);
J3chi2 = zeros(size(delJ3));
J3Iter = zeros(size(delJ3));
for i=1:length(delJ3)
    J1 = fitStr.x(1); J2 = fitStr.x(2);...
        J3 = fitStr.x(3) + delJ3(i);
    
    MnF2Fit = spinw;
    MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
        [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
    MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0],...
        'label', {'MMn2' 'F'}, 'color', {'red' 'yellow'});
    MnF2Fit.gencoupling('maxDistance', intDist);
    MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
    MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
    MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
    MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
    MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
    MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
    MnF2Fit.coupling.rdip = intDist; % Dipole-dipole interaction
    MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);
    
    x1 = [J1 J2];
    MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2'}, 'init', true)
    
    MnF2Fit.table('matrix')
    
    pref = swpref;
    pref.tid = 0;
    
    par_fit          = struct;
    par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';
    
    par_fit.Evect    = linspace(0, 10, 501);
    par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
                        'mat', {'J1' 'J2'}, 'init', true);
    par_fit.xmin      = [-0.6 -0.6];
    par_fit.xmax      = [0.6 0.6];
    par_fit.x0        = x1;
    par_fit.plot      = false;
    par_fit.hermit    = false;
    par_fit.optimizer = 'pso';
    par_fit.maxiter   = 1e2;
    par_fit.nrun      = 1;
    par_fit.TolX      = 10^(-4);
    par_fit.TolFun    = 10^(-5);
    
    fitStrIter = MnF2Fit.fitspec(par_fit);

    J3Iter(i) = fitStr.x(3) + delJ3(i);
    J3chi2(i) = fitStrIter.redX2;
    
    disp(['i: ' num2str(i)])
    disp(['p: ' num2str(fitStrIter.x)])
    disp(['chi2R: ' num2str(fitStrIter.redX2)])
end
chi2RThr = (1 + 1/dof)*min(J3chi2); % Threshold reduced chi-squared
indClose = find(J3chi2<chi2RThr); % J3's within the threshold
sigJ3 = (J3Iter(indClose(end)+1) - J3Iter(indClose(1)-1))/2; % Go to points just beyond intercept with threshold, subtract, divide by two to slightly overestimate uncertainty
figure
hold on
scatter(J3Iter, J3chi2)
yline(chi2RThr)
xlabel('J3 (meV)')
ylabel('Chi2R')
title(['J3: ' num2str(fitStr.x(3)) '. sigJ3: ' num2str(sigJ3)])
box on
hold off
saveas(gcf, strcat('U:/Data/Figures/MnF2/Fitting/UncJ3Dip.pdf'))
disp(['J3: ' num2str(fitStr.x(3))])
disp(['sigJ3: ' num2str(sigJ3)])

%% Global Fit: Construct the best-fit Hamiltonian

J1 = fitStr.x(1); J2 = fitStr.x(2); J3 = fitStr.x(3);

MnF2Fit = spinw;
MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
    [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label',...
    {'MMn2' 'F'}, 'color', {'red' 'yellow'});
MnF2Fit.gencoupling('maxDistance', intDist);
MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
MnF2Fit.addmatrix('label', 'J3', 'value', J3, 'color', 'r');
MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
MnF2Fit.addcoupling('mat', 'J3', 'bond', 3);
MnF2Fit.coupling.rdip = intDist; % Dipole-dipole interaction
MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);

%% Global Fit 2: Read in the data

% Plot over the experimental data after loading the data.
T = sw_readtable('U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt');

% List of Q values from data.
Q = [[T(:).QH]; [T(:).QK]; [T(:).QL]]';
% List of energy values from data.
w = [T(:).EN1];
werr = [T(:).s1];
% List of intensity values from data.
inten = [T(:).I1];
% Indices corresponding to each slice
% HHL
ind1 = 1:5;
ind2 = 6:25;
ind3 = 26:39;
ind4 = 40:55;
ind5 = 56:66;
% H0L
ind6 = 67:83;
ind7 = 84:100;
ind8 = 101:111;
ind9 = 112:121;
% HK0
ind10 = 122:141;

save('U:\Data\MATLAB\MnF2\HamiltonianFitDip.mat')

%% Figure 2: Plot the simulated high-symmetry paths (tutorial 35)

spec = MnF2Fit.spinwave({[0,0,1] [0.5,0.5,1] [0.5,0.5,0.5] [1,1,0.5] [1,1,0] 223}, 'formfact', true,... 
    'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(0, 9, 1000), 'component',...
    'SPerp', 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmE);
figure
sw_plotspec(spec, 'mode', 4, 'qlabel', {'\Gamma' 'M' 'A' 'Z' '\Gamma'})
ylim([0 8.4]);
clim([0 5.5]);
colormap(gca, 'turbo');
legend off;
hold on
% Gamma to M data
wind1 = [5 57:65]; % Indices from datafile
Qabs1 = sqrt(sum((bsxfun(@minus,Q(wind1,:),[0 0 1])*MnF2Fit.rl).^2,2)); % Distance starting from 0
dist1 = sqrt(sum(([0.5 0.5 0]*MnF2Fit.rl).^2)); % Total distance by end
% M to A data
wind2 = flip(16:26);
Qabs2 = sqrt(sum((bsxfun(@minus,Q(wind2,:),[0.5 0.5 1])*MnF2Fit.rl).^2,2)) + dist1;
dist2 = dist1 + sqrt(sum(([0 0 0.5]*MnF2Fit.rl).^2));
% A to Z data
wind3 = 48:56;
Qabs3 = sqrt(sum((bsxfun(@minus,Q(wind3,:),[0.5 0.5 0.5])*MnF2Fit.rl).^2,2)) + dist2;
dist3 = dist2 + sqrt(sum(([0.5 0.5 0]*MnF2Fit.rl).^2));
% Z to Gamma data
wind4 = flip(27:35);
Qabs4 = sqrt(sum((bsxfun(@minus,Q(wind4,:),[1 1 0.5])*MnF2Fit.rl).^2,2)) + dist3;
wabs = w([wind1, wind2, wind3, wind4])';
wabserr = werr([wind1, wind2, wind3, wind4])';
Qabs = [Qabs1; Qabs2; Qabs3; Qabs4];
intAbs = inten([wind1, wind2, wind3, wind4])';
for i=1:length(wabs) % Loop through to vary marker size based on intensity
    plot3(Qabs(i), wabs(i), 1e2+wabs(i)*0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 3/max(intAbs)*intAbs(i)+3)
end
plot3([Qabs, Qabs]', [-wabserr, wabserr]'+wabs', [1e2+wabs*0, 1e2+wabs*0]', 'w-', 'MarkerFaceColor', 'w', 'LineWidth', 1)
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 10, 5]);
for i = 1:length(ax.XTick) % Add vertical lines at each Q that are solid white lines
    xline(ax.XTick(i), 'Color', 'w', 'LineWidth', 1.5)
end
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/HighSymHHL.pdf'))

spec = MnF2Fit.spinwave({[1,0,0] [1.5,0,0] [1.5,0,0.5] [1,0,0.5] 223}, 'formfact', true,... 
    'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(0, 9, 1000), 'component',...
    'SPerp', 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmE);
figure
sw_plotspec(spec, 'mode', 4, 'qlabel', {'\Gamma' 'X' 'R' 'Z'})
ylim([0 8.4]);
clim([0 4]);
colormap(gca, 'turbo');
legend off;
hold on
% Gamma to X data
wind1 = 74:84; % Indices from datafile
Qabs1 = sqrt(sum((bsxfun(@minus,Q(wind1,:),[1 0 0])*MnF2Fit.rl).^2,2)); % Distance starting from 0
dist1 = sqrt(sum(([0.5 0 0]*MnF2Fit.rl).^2)); % Total distance by end
% X to R data
wind2 = 113:122;
Qabs2 = sqrt(sum((bsxfun(@minus,Q(wind2,:),[1.5 0 0])*MnF2Fit.rl).^2,2)) + dist1;
dist2 = dist1 + sqrt(sum(([0 0 0.5]*MnF2Fit.rl).^2));
% R to Z data
wind3 = flip(91:101);
Qabs3 = sqrt(sum((bsxfun(@minus,Q(wind3,:),[1.5 0 0.5])*MnF2Fit.rl).^2,2)) + dist2;
wabs = w([wind1, wind2, wind3])';
wabserr = werr([wind1, wind2, wind3])';
Qabs = [Qabs1; Qabs2; Qabs3];
intAbs = inten([wind1, wind2, wind3])';
for i=1:length(wabs) % Loop through to vary marker size based on intensity
    plot3(Qabs(i), wabs(i), 1e2+wabs(i)*0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 3/max(intAbs)*intAbs(i)+3)
end
plot3([Qabs, Qabs]', [-wabserr, wabserr]'+wabs', [1e2+wabs*0, 1e2+wabs*0]', 'w-', 'MarkerFaceColor', 'w', 'LineWidth', 1)
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 7.5, 5]);
set(gca, 'YTickLabel', [])
for i = 1:length(ax.XTick) % Add vertical lines at each Q that are solid white lines
    xline(ax.XTick(i), 'Color', 'w', 'LineWidth', 1.5)
end
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/HighSymH0L.pdf'))

spec = MnF2Fit.spinwave({[1,0.5,0] [1.5,0.5,0] 223}, 'formfact', true,... 
    'hermit', true);
spec = sw_egrid(spec, 'Evect', linspace(0, 9, 1000), 'component',...
    'SPerp', 'T', 1.6, 'ImagChk', true);
spec = sw_instrument(spec, 'dE', fwhmE);
figure
sw_plotspec(spec, 'mode', 4, 'qlabel', {'X' 'M'})
ylim([0 8.75]);
clim([0 1]);
colormap(gca, 'turbo');
legend off;
hold on
% X to M data. One path so plots against H rather than Qabs.
wind = 132:142;
intAbs = inten(wind)';
for i=1:length(wind) % Loop through to vary marker size based on intensity
    plot3(Q(wind(i),1)-1, w(wind(i)), 1e2+w(wind(i))*0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 3/max(intAbs)*intAbs(i)+3)
end
plot3([Q(wind,1)-1, Q(wind,1)-1]', [-werr(wind); werr(wind)]+w(wind), [1e2+w(wind)*0; 1e2+w(wind)*0], 'w-', 'MarkerFaceColor', 'w', 'LineWidth', 1)
xticks([0 0.5])
xticklabels({'X' 'M'}) % qlabel fails for two points
title('')
ax = gca;
ax.FontSize = 16;
xlabel('')
ylabel('')
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 1, 2.6, 5]);
ylabel('\rm\Delta\it{E}\rm{ (meV)}')
for i = 1:length(ax.XTick) % Add vertical lines at each Q that are solid white lines
    xline(ax.XTick(i), 'Color', 'w', 'LineWidth', 1.5)
end
disp(['Max delta E: ', num2str(max(spec.omega(1,:) - spec.omega(2,:)))])
saveas(gcf, strcat('U:/Data/Figures/MnF2/Manuscript/Fig2/HighSymHK0.pdf'))

%% Use SpinW to calculate analytical expression for dispersion. Tutorial 17

clear, close all

% J1 = 2*-0.028; J2 = 2*0.152; J3 = 2*0.004; Dc = -0.0; % Nikotin 1969
% J1 = -0.06668; J2 = 0.30211; J3 = -0.0041079; Dc = -0.027351;% Previous fitting attempt
% J1 = 2*-0.028; J2 = 2*0.152; J3 = 2*0.000; Dc = -0.091; % Okazaki 1963
J1 = 2*-0.0354; J2 = 2*0.1499; J3 = 2*0.000; Dc = -0.131; % CAMEA 2023
% J1 = 2*-0.094; J2 = 2*0.126; J3 = 2*0.000; Dc = -0.115; % Test

MnF2Fit = spinw;
MnF2Fit.symbolic(true)
MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
    [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label',...
    {'MMn2' 'F'}, 'color', {'red' 'yellow'});
MnF2Fit.gencoupling('maxDistance', 10);
% MnF2Fit.table('bond',1:10)
MnF2Fit.addmatrix('label', 'J1', 'value', -sym('J1', ["positive" "real"]), 'color', 'b');
MnF2Fit.addmatrix('label', 'J2', 'value', sym('J2', ["positive" "real"]), 'color', 'g');
D_mat = @(D) D*diag([0 0 1]);
MnF2Fit.addmatrix('label', 'D', 'value', -sym(D_mat, ["positive" "real"]))
disp('Symbolic matrix value from symbolic input:')
MnF2Fit.matrix.mat
MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
MnF2Fit.addaniso('D');
MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -1 1]);
disp('Ground state energy meV/spin:')
MnF2Fit.energy
spec = MnF2Fit.spinwave();
pretty(spec.omega)

%% Fit to the analytical dispersion relation as given in Yamani 2010. For now J1, J2 though can do J3 from Okazaki 1964 expression.

clear, close all

% J1 = -0.028; J2 = 0.152; Dc = -0.0; % Nikotin 1969
% J1 = -0.028; J2 = 0.152; Dc = -0.091; % Okazaki 1963
J1 = -0.0354; J2 = 0.1499; Dc = -0.131; % CAMEA 2023
% J1 = -0.094; J2 = 0.126; Dc = -0.115; % Test
% J1 = -0.06668; J2 = 0.30211; Dc = -0.027351; % Previous fitting attempt

% Instrumental E-resolution extracted from MJOLNIR for HHL along (001).
% Will be fit to polynomial in sw_instrument.
sigE = [0.3090, 0.0947; 1.0043, mean([0.0917, 0.0948]);...
    2.0002, mean([0.0889, 0.0919]); 2.9962, mean([0.1340, 0.1372]);...
    3.9921, mean([0.1321, 0.1351]); 5.0069, mean([0.1793, 0.1827]);...
    6.0028, mean([0.1779, 0.1811]); 6.9988, mean([0.2276, 0.2310]);...
    7.9947, mean([0.2264, 0.2298])];
fwhmE = sigE; fwhmE(:,2) = sigE(:,2)*2*sqrt(2*log(2));

MnF2Fit = spinw;
MnF2Fit.genlattice('lat_const', [4.8736 4.8736 3.2998], 'angled',...
    [90 90 90], 'sym', 136); % Jauch 1988 11 K X-ray
MnF2Fit.addatom('r', [0 0.30443; 0 0.30443; 0 0], 'S', [5/2 0], 'label',...
    {'MMn2' 'F'}, 'color', {'red' 'yellow'});

MnF2Fit.gencoupling('maxDistance', 10);
% MnF2Fit.table('bond',1:10)

MnF2Fit.addmatrix('label', 'J1', 'value', J1, 'color', 'b');
MnF2Fit.addmatrix('label', 'J2', 'value', J2, 'color', 'g');
MnF2Fit.addmatrix('label', 'Dc', 'value', diag([0 0 Dc]))

MnF2Fit.addcoupling('mat', 'J1', 'bond', 1);
MnF2Fit.addcoupling('mat', 'J2', 'bond', 2);
MnF2Fit.addaniso('Dc');

MnF2Fit.genmagstr('mode', 'direct', 'k', [0 0 0], 'S', [0 0; 0 0; -4.6 4.6]);

% Relax magnetic structure
%MnF2Fit.optmagsteep('nRun', 1e4)
% MnF2Fit.table('mag')

% plot(MnF2Fit, 'range', [1 1 1], 'bondMode', 'line', 'bondLinewidth0', 3)

% Global Fit 1: Prepare fitting to the dispersion (tutorial 35)

x1 = [J1 J2 Dc];
MnF2Fit.matparser('param', x1, 'mat', {'J1' 'J2' 'Dc(3,3)'}, 'init', true)

MnF2Fit.table('matrix')

pref = swpref;
pref.tid = 0;

par_fit          = struct;
par_fit.datapath = 'U:\Data\MATLAB\MnF2\1p5K_0T_FitDisp.txt';

par_fit.Evect    = linspace(0, 10, 501);
% par_fit.func     = @(obj, p) matparser(obj, 'param', p,...
%                     'mat', {'J1' 'J2' 'Dc(3,3)'}, 'init', true);
% zeta = @(obj, p) p(3) + 2*(2*5/2*2*p(1))*sin(qc.*3.2998/2).^2;
% gamma = @() ;
% par_fit.func     = @(obj, p) 
par_fit.xmin      = [-0.6 -0.6 -0.5];
par_fit.xmax      = [0.6 0.6 0.0];
par_fit.x0        = x1;
par_fit.plot      = true;
par_fit.hermit    = false;
par_fit.optimizer = 'pso';
par_fit.maxiter   = 5e2;
par_fit.nrun      = 1;
par_fit.TolX      = 10^(-4);
par_fit.TolFun    = 10^(-5);

% Global Fit 1: Fit to the dispersion (tutorial 35)

fitStr = MnF2Fit.fitspec(par_fit);

disp(['p: ' num2str(fitStr.x)])
disp(['chi2: ' num2str(fitStr.redX2)])
