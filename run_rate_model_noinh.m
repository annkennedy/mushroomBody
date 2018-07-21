realonly = 1;   % realonly == 1 only models glomeruli from which H&C recorded (n = 23),
                % realonly == 0 makes fake data for the remainder of glomeruli (not recommended)

ORN = {}; %olfactory receptor neurons
% odorTotal   = 110;
ORN         = buildORN(ORN, realonly);

% new code to add mixtures of odors ---------------------------------------
% mixinds     = 5:5:110;                  %select odors to make mixtures of
% mixlist     = ORN.odornames(mixinds);
% ORN.rates   = ORN.rates(:,mixinds);     %remove the remainder of the original dataset
% ORN.odornames = mixlist;
% ORN.odorclass = ORN.odorclass(mixinds);
% ORN         = addmixtures(ORN, mixlist);
% -------------------------------------------------------------------------

nOdor = size(ORN.rates,2);

LN={};  %lateral neurons
LN          = buildLN(LN);


for rep = 1
    %set up PN, KC models
    [PN KC] = initialize_cells(ORN,realonly);
    KC.tau_m    = 0.01;
    KC.tau_s    = 0.01;


    %this measures the amount of spontaneous input to KCs
    [~,PN_t,~]  = getPNdynamics(ORN,PN,LN,1);
    PN_t = kron(PN_t,ones(PN.PN_per_glom,1));
    PN_spont    = mean(PN_t(:,(PN.tOn-PN.tmin)/PN.dt/2:(PN.tOn-PN.tmin)/PN.dt),2);
    KC.spont    = KC.wPNKC'*PN_spont;

    %this tunes KC representation sparsity
    clear PN_t;
% fit sparsity with inhibition off
    sp_target = 0.1;
        KC.wInhKC = zeros(KC.ncells,1);
        KC.thr = 1e5;

        count = 0; odorset = 2:3:110;
        KCpks   = zeros(KC.ncells,length(odorset));
        for odorid = odorset
            count = count+1;
            [~,PN_t,~] = getPNdynamics(ORN,PN,LN,odorid);
            PN_stim = kron(PN_t,ones(PN.PN_per_glom,1));

            [KC_rasters,KC_t,~] = getKCratedynamics(PN_stim,KC);
            KCpks(:,count) = max(KC_t')' - KC.spont*2;
        end

        pksort = sort(KCpks(:),'descend');
        thr = pksort(min(ceil(sp_target*KC.ncells*length(odorset)),length(pksort)))
        KC.thr = thr + KC.spont*2;
%


    odorspercell = zeros(KC.ncells,1);
    fr_active   = zeros(nOdor,1);
    KCmean_st   = zeros(KC.ncells,nOdor);
    inh         = [];

    for odorid = 111:nOdor
        [~,PN_t,~]  = getPNdynamics(ORN,PN,LN,odorid);
        PN_t = kron(PN_t,ones(PN.PN_per_glom,1));

        [KC_rasters,KC_t,inh(odorid,:)] = getKCratedynamics(PN_t,KC);

        KCmean = sum(KC_rasters,2);
        fr_active(odorid) = length(nonzeros(KCmean>0))/KC.ncells;

        odorspercell(KCmean>0) = odorspercell(KCmean>0) + 1;

        KCmean_st(:,odorid) = KCmean;
        disp([num2str(odorid) '. ' ORN.odornames{odorid} ': ' num2str(fr_active(odorid)*100) '%']);
    end

    %saves output for use with other functions (eg evan_compare_distances)
    save(['KCsim_withconc' num2str(rep) '.mat'],'PN','KC','fr_active','KCmean_st','inh');
end
