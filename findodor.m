function odorid = findodor(ORN,odorname,varargin)
% find odor index by its name

odorid = find(1-cellfun(@isempty,strfind(ORN.odornames,odorname)));
if(~isempty(varargin))
    if(strcmp(varargin,'exact'))
        odorid = find(strcmp(ORN.odornames,odorname));
    end
end

if(isempty(odorid))
    disp(['warning: could not find odor ' odorname]);
end

if(length(odorid)>1)
    disp(['warning: multiple entries found for ' odorname]);
end