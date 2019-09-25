
% options for the simulation. See get_MB_default_settings for definitions
defopts = get_MB_default_settings();
runOdors = 1:110;
 
% unused code to simulate responses to odor mixtures: 
% defopts.addMixtures = 1;
% ORN = buildORNs({}, opts.realonly);
% [public, private] = find_public_and_private_odors(ORN,3);
% defopts.mixinds     = [public;private];
% defopts.ratios      = 0.2:0.2:0.8;
% runOdors = 1:(110 + (length(defopts.mixinds)*(length(defopts.mixinds)-1)/2) * length(defopts.ratios));

% design the experiments to run:
experiment(1) = defopts;
experiment(2) = defopts;  experiment(2).homeostatic=1;
experiment(3) = defopts;  experiment(3).shuffleORNs=1;
experiment(4) = defopts;  experiment(4).uniformPNKC=1;
experiment(5) = defopts;  experiment(5).useStaticPN=1;
defopts.useAPL=0;
experiment(6) = defopts;
experiment(7) = defopts;  experiment(7).homeostatic=1;
experiment(8) = defopts;  experiment(8).shuffleORNs=1;
experiment(9) = defopts;  experiment(9).uniformPNKC=1;
experiment(10) = defopts; experiment(10).useStaticPN=1;

% and run them!
for i=1:length(experiment)
    fprintf('Running experiment %d/%d...\n', i,length(experiment));
    opts = experiment(i) %load the option settings for this experiment

    
    for rep = 1:opts.nreps
        fprintf('rep %d/%d\n',rep,opts.nreps);
        experiment(i).results{rep} = run_rate_model(runOdors,opts);
    end
end

