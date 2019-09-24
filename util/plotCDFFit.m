function h=plotCDFFit(vals)

valFit = [];
for i=unique(vals)
    valFit(end+1) = mean(vals>=i);
end
h=plot(unique(vals),valFit);
