function LN = buildLN(LN)

% time constants (GA/GB fit to RI Wilson, G Laurent (2005) J Neuro)
LN.taum   = 0.01;
LN.tauGA  = 0.1;
LN.tauGB  = 0.4;

% two parameters for divisive inhibition:
LN.inhsc  = 500;
LN.inhadd = 200;

% will refine this in spiking model (will also need to add a GABA depletion
% time constant)
LN.thr    = 1;