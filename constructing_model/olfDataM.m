
%Original H-C data
[xDat,odorclass,odornames,glomnames,glomnames_sophie] = load_HC_data();
nOdor = length(odornames) - 1;

ORND = xDat(1:nOdor, :);

% adds in spontaneous rates:
ORNBack = xDat(nOdor+1, :);
ORND = ORND + ones(nOdor, 1)*ORNBack;

% remove dm3 for some reason
ORND(:, 8) = [];
% ORN(:, 12) = [];
% ORN(:, 14) = [];
% ORN(:, 20) = [];
ORNBack(:, 8) = [];
% ORNBack(:, 12) = [];
% ORNBack(:, 14) = [];
% ORNBack(:, 20) = [];
glomnames(8)=[];