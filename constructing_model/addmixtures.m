function ORN = addmixtures(ORN, mixlist, ratios)
% first find the strings in mixlist among the entries of ORN (exact match
% required), then append the mixtures to the end of the appropriate ORN
% fields
% ratios gives the fraction of odor A. So eg [.2 .4 .6 .8] would be .2:.8,
% .4:.6, .6:.4, and .8:.2.


stiminds = zeros(size(mixlist));

for i = 1:length(mixlist)
    temp = strcmp(ORN.odornames,mixlist{i});
    if(sum(temp)>0)
        stiminds(i) = find(temp);
    end
end

% get all possible mixtures of the listed odors:
[a,b] = meshgrid(stiminds,stiminds);
pairs = [a(:) b(:)];
pairs(pairs(:,1)>=pairs(:,2),:) = [];

% get the ORN responses and new odor names:
newodors = zeros(size(pairs,1),51);
newnames = cell(1,size(pairs,1));
nrat = length(ratios);
for i=1:size(pairs,1)
    for j=1:nrat
        newodors((i-1)*nrat+j,:)   = (ORN.rates(:,pairs(i,1))*ratios(j) + ORN.rates(:,pairs(i,2))*(1-ratios(j)));
        newnames{(i-1)*nrat+j}     = ['mix ' ORN.odornames{pairs(i,1)} num2str(ratios(j)) ...
                                        ' + ' ORN.odornames{pairs(i,2)} num2str(1-ratios(j))];
    end
end

ORN.rates = [ORN.rates newodors'];
ORN.odornames = [ORN.odornames' newnames]';
ORN.odorclass = [ORN.odorclass 13*ones(1,length(newnames))]; %add a 13th odor class for mixtures