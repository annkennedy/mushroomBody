realonly = 0;   % realonly == 1 only models glomeruli from which H&C recorded (n = 23),
                % realonly == 0 makes fake data for the remainder of glomeruli (not recommended)

ORN = {}; %olfactory receptor neurons
ORN         = buildORN(ORN, realonly,5000);

% new code to add mixtures of odors ---------------------------------------
% mixinds     = randperm(110,20);                  %select odors to make mixtures of
% ratios      = 0.1:0.1:0.9;                       %set mixing ratios
% mixlist     = ORN.odornames(mixinds);
% [ORN,pairs] = addmixtures(ORN, mixlist, ratios);
% -------------------------------------------------------------------------

nOdor = size(ORN.rates,2);
odorlist = 1:nOdor;

LN={};  %lateral neurons (inhbition between glomeruli)
LN          = buildLN(LN);


for rep = 1
    %set up PN and KC models
    [PN, KC] = initialize_cells(ORN,realonly);
    KC.tau_m    = 0.01;
    KC.tau_s    = 0.01;

    %first, measure the amount of spontaneous input to KCs
    [~,PN_t,~]  = getPNdynamics(ORN,PN,LN,1);
    PN_t = kron(PN_t,ones(PN.PN_per_glom,1));
    PN_spont    = mean(PN_t(:,(PN.tOn-PN.tmin)/PN.dt/2:(PN.tOn-PN.tmin)/PN.dt),2);
    KC.spont    = KC.wPNKC'*PN_spont;
    disp('Model initialized');

    %next, tune inhibition to set KC representation sparsity
    disp('Fitting APL inhibition (this step is slow)...');
    clear PN_t;
    odorset = 2:3:110;  % code takes a while to run so I don't use all odors
    sp_target = 0.1;    % target sparsity of model KC population
    for i = odorset
        [~,PN_t{i},~] = getPNdynamics(ORN,PN,LN,i);
    end
    KC = fitSparseness(PN,PN_t,odorset,KC,sp_target);
    disp('APL inhibition fit');


    %now, simulate responses to odors specified by odorlist
    odorspercell = zeros(KC.ncells,1);
    fr_active    = zeros(nOdor,1);
    KCmean_st    = zeros(KC.ncells,nOdor);

    for odorid = odorlist
        [~,PN_t,~]  = getPNdynamics(ORN,PN,LN,odorid);
        PN_t = kron(PN_t,ones(PN.PN_per_glom,1));

        [KC_rasters,KC_t,~] = getKCratedynamics(PN_t,KC);

        KCmean = sum(KC_rasters,2);
        fr_active(odorid) = length(nonzeros(KCmean>0))/KC.ncells;

        odorspercell(KCmean>0) = odorspercell(KCmean>0) + 1;

        KCmean_st(:,odorid) = KCmean;
        disp([num2str(odorid) '. ' ORN.odornames{odorid} ': ' num2str(fr_active(odorid)*100) '%']);
    end

    %save output for use with other functions (most of these I haven't included in this code)
%     save(['KCsim_' num2str(rep) '.mat'],'PN','KC','fr_active','KCmean_st','ORN');
end
%

s=svd(KCmean_st);
dimensionality = sum(s)^2/sum(s.^2)
