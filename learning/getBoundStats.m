function Stats = getBoundStats(rON,rONpre, useodors,distractodors,classes)

%returns some set of general stats about the quality of separation between
%trained and untrained odors. Stats is a struct.

genodors  = setdiff(1:length(rON),[useodors distractodors]);
trained   = rON(useodors);
untrained = rON(genodors);

%always define stats relative to max(trained)
Stats.mean_distance = mean(untrained)   - max(trained);
Stats.min_distance  = min(untrained)    - max(trained);
Stats.med_distance  = median(untrained) - max(trained);
Stats.var_unt       = std(untrained);
Stats.change_rON    = mean(abs(rONpre(genodors) - rON(genodors)));

Stats.trainodor_cats = classes(useodors);

if(Stats.min_distance<0)
    Stats.false_pos_count   = mean(untrained<max(trained));
    false_pos = find(untrained<max(trained));
    for c=1:length(unique(classes))
        Stats.false_pos_cats(c) = sum(classes(false_pos)==c);
    end
else
    Stats.false_pos_count = 0;
    Stats.false_pos_cats  = zeros(size(unique(classes)));
end