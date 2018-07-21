function [odorresp,odorclass,odornames,glomnames,glomnames_sophie] = load_HC_data()

pth = fileparts(mfilename('fullpath'));
[~,~,allText] = xlsread([pth '\HC_data_raw.csv']);

glomnames           = allText(1,3:end);
glomnames_sophie    = allText(2,3:end);

odorclass   = [allText{3:end,1}];
odornames   = allText(3:end,2);
odorresp    = cell2mat(allText(3:end,3:end));