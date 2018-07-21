function [PN, KC] = initialize_cells(ORN,realonly)
PN={};  %projection neurons
KC={};  %kenyon cells

%% set up PNs

PN.ncells   = ORN.ncells;
PN.PN_per_glom  = 5;  % number of PNs projecting from each glom. have same mean signal, independent spiking noise/targets
PN.HCList       = ORN.HCList;

% parameters for new PN model, fit via least squares:
PN.Rmax     = 165;
PN.offset   = 2.9410;
PN.tanhsc   = 5.3395;
PN.inhsc    = 368.6631;
PN.inhadd   = 31.4088;

% parameters for dynamic simulation
PN.taum     = 0.01;
PN.Rmax_dyn = 200;          %peak firing rate in dynamic model
PN.dt       = 0.5*10^-3;    %all times are in seconds
PN.tmin     = -0.5;         %when to start simulating (relative to odor onset)
PN.tOn      = 0;            %odor onset
PN.tOff     = 0.5;          %odor offset
PN.tmax     = 0.75;         %when to stop simulating
PN.ref      = 0.003/PN.dt;  %refractory period (fit to mean-variance relationship of spikes per trial in RW data)

PN          = PNsim(PN,ORN,realonly);
PN.rates_sp = zeros(PN.ncells*PN.PN_per_glom,size(PN.rates,2));

%% set up KCs

KC.ncells   = 2000;
KC.sparse   = 0.05;
KC.gain     = 100;

%parameters for dynamics simulation
KC.dt       = 0.5*10^-3;
KC.tmin     = PN.tmin;
KC.tmax     = PN.tmax;
KC.tOn      = PN.tOn;
KC.tOff     = PN.tOff;
KC.taum     = 0.01;

KC.wPNKC    = buildConnects(KC,PN, realonly);
PN.wPNInh   = ones(1,PN.ncells);
KC.wKCInh   = ones(1,KC.ncells)/KC.ncells;
KC.wInhKC   = ones(KC.ncells,1)/10;
KC.divisiveInh = 1; %1 for divisive APL inhibition, 0 for subtractive

% [KC inh] = KCsim(KC,PN);

% dofigs;
