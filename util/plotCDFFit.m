function h=plotCDFFit(vals)

valFit = [];
for i=unique(vals)
    valFit(end+1) = mean(vals<=i);% + mean(vals==i)/2;
end
h=plot(unique(vals),valFit);