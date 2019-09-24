function opts = get_MB_default_settings()

% default options for the simulation in run_rate_model.m:

% realonly==1-> only model glomeruli from H&C dataset (n = 23)
opts.realonly    = 1;

% whether to include APL inhibition
opts.useAPL      = 1;

% whether to include homeostatic threshold tuning
opts.homeostatic = 0;

% decorrelate odors by shuffling glom. identities at ORN level
opts.shuffleORNs = 0;

% uniformPNKC==1 -> don't use Caron et al connectivity data
opts.uniformPNKC = 0;

% average fraction of KCs responding to an odor
opts.sparsity    = 0.1;

% ==1 -> use the Olsen et al PN model instead of the dynamic model
opts.useStaticPN = 0;

% make mixtures of odors (must also provide indicies of odors to mix, + mixing ratios,
% as opts.mixinds and opts.ratios)
opts.addMixtures = 0;

% number of repeats to run with given parameter settings
opts.nreps       = 1;

% show progress during simulation
opts.verbose     = 1;
