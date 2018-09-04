function KC = scalePNKCwts(PN, PN_t, odorset, KC)

count = 0;
PNpks   = zeros(PN.ncells*PN.PN_per_glom,length(odorset));
for odorid = odorset
    count = count+1;
    PN_stim = kron(PN_t{odorid},ones(PN.PN_per_glom,1));

    PNpks(:,count) = max(PN_stim')';
end

KC.wPNKC = bsxfun(@times,KC.wPNKC,100./(100+mean(PNpks')'));
