function [KC_spikes,fr_active,odorspercell] = run_rate_model(odorlist, opts)

    if(~exist('opts','var'))
        opts = get_MB_default_settings();
    end
    
    % build all the model neuron populations
    ORN         = buildORNs({}, opts);
    LN          = buildLNs({});
    PN          = buildPNs(ORN,opts);
    KC          = buildKCs(ORN,PN,opts);
    nOdor       = size(ORN.rates,2);

    %first, measure the amount of spontaneous input to KCs
    [~,PN_t,~]  = getPNdynamics(ORN,PN,LN,1);
    PN_t = kron(PN_t,ones(PN.PN_per_glom,1));
    PN_spont    = mean(PN_t(:,(PN.tOn-PN.tmin)/PN.dt/2:(PN.tOn-PN.tmin)/PN.dt),2);
    KC.spont    = KC.wPNKC'*PN_spont;
    if(opts.verbose)
        disp('Model initialized');
    end
    
    %next, tune inhibition to set KC representation sparsity
    clear PN_sample PN_RW_sample;
    odorset = 2:3:110;     % odors to use for threshold-setting.
    sp_target = opts.sparsity; % target sparsity of model KC population
    for i = odorset
        [~,PN_sample{i},~,PN_RW_sample{i}] = getPNdynamics(ORN,PN,LN,i);
    end
    if(opts.useStaticPN)
        PN_sample = PN_RW_sample; % to test with inputs from the static model
    end
    if(opts.useAPL)
        if(opts.verbose)
            disp('Fitting APL inhibition...');
        end
        if(opts.homeostatic)
            KC = setKCthreshold_bycell(PN,PN_sample, odorset, KC,sp_target*2);
        else
            KC = setKCthreshold(PN,PN_sample, odorset, KC,sp_target*2);
        end
        KC = fitSparseness(PN,PN_sample,odorset,KC,sp_target);
    else
        KC.wInhKC = zeros(KC.ncells,1);
        KC.wKCInh = zeros(1,KC.ncells);
        KC = scalePNKCwts(PN, PN_sample, odorset, KC);
        if(opts.homeostatic)
            KC = setKCthreshold_bycell(PN,PN_sample, odorset, KC,sp_target);
        else
            KC = setKCthreshold(PN,PN_sample,odorset,KC,sp_target);
        end
    end
    if(opts.verbose)
        disp('APL inhibition fit');
    end


    %now, simulate responses to odors specified by odorlist
    odorspercell = zeros(KC.ncells,1);
    fr_active    = zeros(nOdor,1);
    KC_spikes    = zeros(KC.ncells,nOdor);
    for odorid = odorlist
		if(opts.useFiringRateModel)
			[~,PN_t,~,PN_RW]  = getPNdynamics(ORN,PN,LN,odorid);
			if(opts.useStaticPN)
				PN_t = PN_RW;
			end
			PN_t = kron(PN_t,ones(PN.PN_per_glom,1));
		else
			PN_t = getPNspikes(PN,PN_t,ntrials);
		end

        [KC_rasters,~,~] = getKCdynamics(PN_t,KC);
        KCmean                  = sum(KC_rasters,2);
        fr_active(odorid)       = length(nonzeros(KCmean>0))/KC.ncells;
        odorspercell(KCmean>0)  = odorspercell(KCmean>0) + 1;
        KC_spikes(:,odorid)     = KCmean;
        
        if(mod(odorid,10)==0 && opts.verbose)
            disp([num2str(odorid) '. ' ORN.odornames{odorid} ': ' num2str(fr_active(odorid)*100) '%']);
        end
    end
