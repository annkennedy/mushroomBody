function KC = fitSparseness(PN,PN_t, odorset, KC,sp_target)

KC.wInhKC = zeros(KC.ncells,1);
% KC.thr    = (-1750*sp_target + 410 + 150*(sp_target<0.01)) + KC.spont*2;
KC.thr = 1e5;

count = 0;
KCpks   = zeros(KC.ncells,length(odorset));
for odorid = odorset
    count = count+1
    PN_stim = kron(PN_t{odorid},ones(PN.PN_per_glom,1));

    [KC_rasters,KC_t,~] = getKCratedynamics(PN_stim,KC);
    KCpks(:,count) = max(KC_t')' - KC.spont*2;
end

pksort = sort(KCpks(:),'descend');
thr = pksort(min(ceil(sp_target*2*KC.ncells*length(odorset)),length(pksort)))
KC.thr = thr + KC.spont*2;
mean(mean(bsxfun(@gt,KCpks,thr)))



if(sp_target<1)
    sp = 0; count = 0;
    KC.wInhKC = ones(KC.ncells,1) * (2*ceil(-log(sp_target)));
    KC.wKCInh = (1/KC.ncells) * ones(1,KC.ncells) * (2*ceil(-log(sp_target)));

    while (abs(sp-sp_target)>0.1*sp_target)&(count<15)
        count = count+1

        odorspercell = zeros(KC.ncells,1);
        fr_active    = zeros(length(odorset),1);
        KCmean_st    = zeros(KC.ncells,length(odorset));
        for odorid = odorset
            PN_stim = kron(PN_t{odorid},ones(PN.PN_per_glom,1));

            [KC_rasters,~,~] = getKCratedynamics(PN_stim,KC);

            KCmean = sum(KC_rasters,2);
            fr_active(odorid) = length(nonzeros(KCmean>0))/KC.ncells;
            odorspercell(KCmean>0) = odorspercell(KCmean>0) + 1;
            KCmean_st(:,odorid) = KCmean;
        end

        sp = mean(mean(KCmean_st(:,odorset)>0));
        [sp sp-sp_target]
        KC.wInhKC = KC.wInhKC + (sp-sp_target)*10/sp_target;
        KC.wKCInh = KC.wKCInh + (sp-sp_target)*10/sp_target/KC.ncells;
    end
end













