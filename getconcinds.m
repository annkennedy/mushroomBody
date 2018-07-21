function concinds = getconcinds(ORN,odortype)
% looks for concentration data from the Hallem & Carlson dataset

if(strcmp(odortype,'fruit'))
    temp = ORN.odornames(find(1-cellfun(@isempty,strfind(ORN.odornames,'pure'))));
elseif(strcmp(odortype,'all'))
    temp = ORN.odornames(find(1-cellfun(@isempty,strfind(ORN.odornames,' -2'))));
else
    temp = ORN.odornames(find(1-cellfun(@isempty,strfind(ORN.odornames,'-8'))));
end

concinds = [];
for i=1:length(temp)
    name = temp{i}(1:max(strfind(temp{i},' ')));
    
    add_data = find(1-cellfun(@isempty,strfind(ORN.odornames,name)));
    if(strfind(name,'apple'))
        drop = find(1-cellfun(@isempty,strfind(ORN.odornames,'pineapple')));
        add_data = setdiff(add_data,drop);
    end
    if(isempty(strfind(name,'mix')))
        inds = cellfun(@isempty,strfind(ORN.odornames(add_data),'mix'));
        concinds = [concinds; add_data(inds)];
    end
    
end
concinds = concinds';
concinds = concinds(:);