function PN = buildPNs(ORN,opts)

PN={};  %projection neurons

PN.ncells   = ORN.ncells;
PN.PN_per_glom  = 5;  % number of PNs projecting from each glom. have same mean signal, independent spiking noise/targets
PN.HCList       = ORN.HCList;

% parameters for new PN model, fit via least squares:
PN.Rmax     = 165;
PN.offset   = 2.9410;
PN.tanhsc   = 5.3395;
PN.inhsc    = 368.6631;
PN.inhadd   = 31.4088;

% parameters for dynamics simulation
PN.tau_m    = 0.01;         % membrane time constant
PN.Rmax_dyn = 200;          % peak firing rate in dynamic model
PN.dt       = 0.5*10^-3;    % simulation timestep in seconds
PN.tmin     = -0.5;         % when to start simulating (relative to odor onset)
PN.tOn      = 0;            % odor onset
PN.tOff     = 0.5;          % odor offset
PN.tmax     = 0.75;         % when to stop simulating
PN.ref      = 0.003/PN.dt;  % refractory period (fit to mean-variance relationship of spikes per trial in RW data)

PN          = PNsim(PN,ORN,opts.realonly); % get PN firing rate dynamics for use with KC
PN.rates_sp = zeros(PN.ncells*PN.PN_per_glom,size(PN.rates,2));
