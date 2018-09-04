function KC = setKCthreshold_bycell(PN,PN_t, odorset, KC,sp_target)

KC.wInhKC = zeros(KC.ncells,1);
KC.thr = 1e5;

count = 0;
KCpks   = zeros(KC.ncells,length(odorset));
for odorid = odorset
    count = count+1;
    PN_stim = kron(PN_t{odorid},ones(PN.PN_per_glom,1));

    [KC_rasters,KC_t,~] = getKCratedynamics(PN_stim,KC);
    KCpks(:,count) = max(KC_t')' - KC.spont*2;
end

% set each cell's threshold to respond to roughly sp_target % of odors
for i = 1:2000
    pksort = sort(KCpks(i,:),'descend');
    thr = pksort(round(sp_target*length(pksort)));
    KC.thr(i,1) = thr + KC.spont(i)*2;
end
