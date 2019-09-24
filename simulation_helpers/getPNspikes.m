function raster = getPNspikes(PN,PN_t,ntrials)

dt      = PN.dt;
time    = PN.tmin+dt:dt:PN.tmax;

raster = cell(ntrials,1);

for i=1:ntrials

    raster{i} = zeros(PN.ncells*PN.PN_per_glom,length(time));
    rate = kron(PN_t*PN.dt/.5,ones(PN.PN_per_glom,1));

    for t=1:length(time)

        win = max(t-PN.ref:t-1,1);
        rate(:,t) = rate(:,t).*(sum(raster{i}(:,win)')'==0); %leave rate unchanged only if there hasn't been a spike lately
        raster{i}(:,t) = poissrnd(rate(:,t));

    end
end