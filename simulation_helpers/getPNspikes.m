function raster = getPNspikes(PN,PN_t,ntrials)

dt      = PN.dt;
time    = PN.tmin+dt:dt:PN.tmax;

raster = cell(ntrials,1);

for i=1:ntrials

    raster{i} = zeros(PN.ncells*PN.PN_per_glom,length(time));
    rate = PN_t*PN.dt/.5;

    for t=1:length(time)

        win = max(t-PN.ref:t-1,1);
        rate(:,t) = rate(:,t).*(sum(raster{i}(:,win)')'==0); %leave rate unchanged if there hasn't been a spike lately
        raster{i}(:,t) = poissrnd(rate(:,t));

    end
end