temp = [(1:110)' sort(ORN.rates(ORN.HCList,1:110),'descend')'];

numodors=10;


temp_srt = sortrows(temp,9);

clc;
disp('public odors:');
public = temp_srt(end-numodors+1:end,1);
ORN.odornames(public)'

temp = sort(ORN.rates(ORN.HCList,1:110),'descend')';

temp = [(1:110)' temp(:,1)-temp(:,4) temp];
temp_srt = sortrows(temp,2);


disp('private odors:')
private = temp_srt(end-numodors:end-1,1);
ORN.odornames(private)'