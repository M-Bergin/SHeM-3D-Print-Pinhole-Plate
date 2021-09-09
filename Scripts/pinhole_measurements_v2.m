% Script to plot experimental and simulated z scans with and without
% shielding around the pinhole plate

% Setup parameters
empty_ = cell([1 2]);
all_ = 1:2;

WD = empty_;
Bright = empty_;
Dark = empty_;
Reference = empty_;

Ie1_model = empty_;

% Modelled scattering efficiency

eta_ss = empty_;
eta_sm = empty_;
eta_es = empty_;
eta_em = empty_;

tau = 1; %Set to 1 to consider the raw count rate

G_ion = 5.7;

%% Load data without Shield
%Experimental data
data = readmatrix('Pinhole Plate Measurements.xlsx', ...
                  'Sheet', '2020-07-07 z Scan', ...
                  'Range', 'A2:F37', ...
                  'UseExcel', false);
              
WD{1} = 2.56 - data(:,1)/1000;
Bright{1} = data(:,5);
Dark{1} = data(:,4);
Reference{1} = data(:,6);

data = readmatrix('Pinhole Plate Measurements.xlsx', ...
                  'Sheet', '2020-07-07 Sensitivity', ...
                  'Range', 'A2:B7', ...
                  'UseExcel', false);
              
Ie1_model{1} = fit(data(:,1)*G_ion*100, data(:,2)/(30*tau), 'poly1');

%Simulated data
sim = load('skeleton (high density).mat');
WD_m = sim.wd;
eta_ss{1} = [sim.supersonic_results.single]./(sim.n_rays - [sim.supersonic_results.killed]);
eta_sm{1} = [sim.supersonic_results.multiple]./(sim.n_rays - [sim.supersonic_results.killed]);
eta_es{1} = [sim.effusive_results.single]./(sim.n_rays - [sim.effusive_results.killed]);
eta_em{1} = [sim.effusive_results.multiple]./(sim.n_rays - [sim.effusive_results.killed]);

%% Load data with Shield
%Experimental data
data = readmatrix('Pinhole Plate Measurements.xlsx', ...
                  'Sheet', '2020-07-01 z Scan II', ...
                  'Range', 'B2:G37', ...
                  'UseExcel', false);
              
% WD{2} = 2.22 - data(:,1)/1000;
WD{2} = 2.59 - data(:,1)/1000;
Bright{2} = data(:,5);
Dark{2} = data(:,4);
Reference{2} = data(:,6);

data = readmatrix('Pinhole Plate Measurements.xlsx', ...
                  'Sheet', '2020-07-01 Sensitivity', ...
                  'Range', 'A2:B7', ...
                  'UseExcel', false);
              
Ie1_model{2} = fit(data(:,1)*G_ion*100, data(:,2)/(30*tau), 'poly1');

%Simulated data
sim = load('shield (high density).mat');
eta_ss{2} = [sim.supersonic_results.single]./(sim.n_rays - [sim.supersonic_results.killed]);
eta_sm{2} = [sim.supersonic_results.multiple]./(sim.n_rays - [sim.supersonic_results.killed]);
eta_es{2} = [sim.effusive_results.single]./(sim.n_rays - [sim.effusive_results.killed]);
eta_em{2} = [sim.effusive_results.multiple]./(sim.n_rays - [sim.effusive_results.killed]);

%% Post-process data
Is = empty_;
Ie = empty_;
for i = all_
    r = Reference{i};
    cal = r(1)./r;
    Ie{i} = 2*cal./Dark{i};
    Is{i} = 2*cal./Bright{i} - Ie{i};
    Ie{i} = Ie{i} / tau;
    Is{i} = Is{i} / tau;
end

%% Scale simulated data
eta_s = arrayfun(@(i)eta_ss{i}+eta_sm{i}, all_, 'UniformOutput', false);
% eta_e = arrayfun(@(i)eta_es{i}+eta_em{i}, all_, 'UniformOutput', false);

f = @(Ns,WD0,x) Ns.*spline(WD_m, eta_s{1}, x-WD0);

[model_s,gof_s] = fit(WD{1}, Is{1}, f, ...
    'StartPoint', [5e5/tau 0], ...
    'Lower', [1e5/tau -.5], ...
    'Upper', [1e6/tau  .5]);

ci_s = confint(model_s);

%% Plot the data
figure
hold on
box on
plot(WD{1}-model_s.WD0, Is{1}, '.', 'MarkerSize', 10,'Color',[150 150 150]/255)
plot(WD_m, model_s.Ns*eta_s{1},'LineWidth',1,'Color',[150 150 150]/255)%[0, 0.4470, 0.7410])
% plot(WD_m, ci(1,1)*eta_s{1}, '--')
% plot(WD_m, ci(2,1)*eta_s{1}, '--')

plot(WD{2}-model_s.WD0, Is{2}, '.', 'MarkerSize', 10,'Color',[229 156 0]/255	)%[0.8500, 0.3250, 0.0980])
plot(WD_m, model_s.Ns*eta_s{2},'LineWidth',1,'Color',	[229 156 0]/255)%[0.8500, 0.3250, 0.0980])

xlim([0 4.2])
set(gca,'FontSize',12,'LineWidth',1)
xlabel('Perpendicular sample distance/mm')
ylabel('Primary beam signal/Hz')

% exportgraphics(gcf,'..\Figures\Shielding_z_scan.eps')
% savefig(['..\Figures\Shielding_z_scan.fig'])
