function low = scaleStd(mean, std, matrix, fig_cut)

low = zeros(length(mean),1)+matrix(fig_cut+10, 4);

for n=1:length(mean)
    if ((mean(n)-std(n)) > 0)
        low(n) = mean(n)-std(n);
    end
end

end