function datasheetclean = trim(datasheet)
% For some reason, MatLab reduplicates the last line of some of the data sheets produced
% by <peakdet.m>. This script trims data sheets to remove this problem.
a = length(nonzeros(datasheet(:,1)));

if a > 1
    if datasheet(a,1) == datasheet(a - 1,1)
        datasheet(a,:) = 0;
        disp('Corrected a data sheet whose last line was identical to the line before last.')
    end
end
datasheetclean = datasheet;