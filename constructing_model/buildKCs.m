function KC = buildKCs(ORN,PN,opts)

KC = {};

KC.ncells   = 2000;
KC.sparse   = 0.05;
KC.gain     = 100;

%parameters for dynamics simulation
KC.tau_m    = 0.01;         % membrane time constant
KC.tau_s    = 0.01;         % synaptic time constant
KC.dt       = 0.5*10^-3;    % smulation timestep in seconds
KC.tmin     = PN.tmin;      % inherited for convenience
KC.tmax     = PN.tmax;
KC.tOn      = PN.tOn;
KC.tOff     = PN.tOff;

KC.wPNKC    = buildConnects(KC,PN,opts.realonly);
if(opts.uniformPNKC) % generate new weight matrix
    populated = ORN.HCList * PN.PN_per_glom;
    populated = sort([populated-4 populated-3 populated-2 populated-1 populated]);
    for i = 1:2000
        KC.wPNKC(populated,i) = KC.wPNKC(populated(randperm(length(populated))),i);
    end
end

PN.wPNInh   = ones(1,PN.ncells);
KC.wKCInh   = ones(1,KC.ncells)/KC.ncells;
KC.wInhKC   = ones(KC.ncells,1)/10;
KC.divisiveInh = 1; %1 for divisive APL inhibition, 0 for subtractive