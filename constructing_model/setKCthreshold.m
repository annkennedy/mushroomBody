function KC = setKCthreshold(PN,PN_t, odorset, KC,sp_target)

KC.wInhKC = zeros(KC.ncells,1);
KC.thr = 1e5;

count = 0;
KCpks   = zeros(KC.ncells,length(odorset));
for odorid = odorset
    count = count+1;
    PN_stim = kron(PN_t{odorid},ones(PN.PN_per_glom,1));

    [KC_rasters,KC_t,~] = getKCdynamics(PN_stim,KC);
    KCpks(:,count) = max(KC_t')' - KC.spont*2;
end

pksort = sort(KCpks(:),'descend');
thr = pksort(min(ceil(sp_target*KC.ncells*length(odorset)),length(pksort)));
KC.thr = thr + KC.spont*2;
