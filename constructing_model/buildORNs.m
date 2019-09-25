function ORN = buildORNs(ORN, opts, varargin)

if(~exist('opts','var'))
	opts = get_MB_default_settings();
end

ORN.HCList  = [7, 17, 46, 12, 8, 20, 5, 39, 6, 45, 21, 29, 33, 22, 15, 24, 40, 34, 23, 48, 16, 28, 49]; %the glomeruli H&C record from

olfDataM;
ORND = ORND - ones(nOdor, 1)*ORNBack; %remove the spontaneous activity

nP = 51;
ORNFull = zeros(nOdor, nP);
ORNFull(:, 1:23) = ORND;

if(~opts.realonly) %make fake data for some of the glomeruli
    for i=24:nP
        for j=1:nOdor
            ORNFull(j, i) = ORND(j, randi(23, 1, 1));
        end
    end
end

%deal ORNFull into the appropriate glomeruli, using H&C's labeling
ORN.rates = zeros(nOdor, nP);
ORN.rates(:,ORN.HCList) = ORNFull(:,1:23);
ORN.rates(:,setdiff(1:51,ORN.HCList)) = ORNFull(:,24:end); %either zeros (realonly==1) or fake data (realonly==0)

ORN.rates       = ORN.rates';
ORN.odornames   = odornames(1:end-1);
ORN.odorclass   = odorclass(1:end-1);
ORN.glomnames   = glomnames;
ORN.glomnames_sophie = glomnames_sophie;
ORN.ncells      = size(ORN.rates,1);

if(~opts.realonly)
    ORN.spont   = ORNBack(randi(23,51,1));
else
    ORN.spont       = zeros(51,1);
end
ORN.spont(ORN.HCList) = ORNBack;

ORN.taum        = 0.01; %membrane time constant for dynamic model

%generate fake odors
if(~isempty(varargin))
    nTotal = varargin{1};
    ORN.rates = ORN.rates(:,1:110);
    for cellnum = 1:length(ORN.HCList)
        ORN.rates(ORN.HCList(cellnum),111:nTotal) = ORN.rates(ORN.HCList(cellnum),randi(110,nTotal-110,1));
    end
    ORN.odorclass(111:nTotal)=max(ORN.odorclass)+1;
    for i=111:nTotal
        ORN.odornames{i} = ['synth ' num2str(i-110)];
    end
end

% optional: make odor mixtures
if(opts.addMixtures)
    mixlist     = ORN.odornames(opts.mixinds);
    ORN = addmixtures(ORN, mixlist, opts.ratios);
end

% optional: shuffles glomerulus identities to decorrelate odors
if(opts.shuffleORNs)
    for i=1:size(ORN.rates,2)
        ORN.rates(ORN.HCList,i) = ORN.rates(ORN.HCList(randperm(23)),i);
    end
end