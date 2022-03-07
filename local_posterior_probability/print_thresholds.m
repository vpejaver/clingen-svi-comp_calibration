function [] = print_thresholds (tp, tb, dtp, dtb)


fprintf(1, '\n\nBenign\n');
for i = 1 : size(tb, 2)
    if i == 1
        fprintf(1, 'VS: ');
    elseif i == 2
        fprintf(1, 'ST: ');
    elseif i == 3
        fprintf(1, 'MO: ');
    else
        fprintf(1, 'SU: ');
    end
        
    if isnan(tb(i))
        fprintf(1, '------; ');
    else
        fprintf(1, '%.4f; ', tb(i));
    end
end

fprintf(1, '\n\nDiscounted Benign\n');
for i = 1 : size(tb, 2)
    if i == 1
        fprintf(1, 'VS: ');
    elseif i == 2
        fprintf(1, 'ST: ');
    elseif i == 3
        fprintf(1, 'MO: ');
    else
        fprintf(1, 'SU: ');
    end
        
    if isnan(dtb(i))
        fprintf(1, '------; ');
    else
        fprintf(1, '%.4f; ', dtb(i));
    end
end

fprintf(1, '\n\nPathogenic\n');
for i = size(tp, 2) : -1 : 1
    if i == 1
        fprintf(1, 'VS: ');
    elseif i == 2
        fprintf(1, 'ST: ');
    elseif i == 3
        fprintf(1, 'MO: ');
    else
        fprintf(1, 'SU: ');
    end
    
    if isnan(tp(i))
        fprintf(1, '------; ');
    else
        fprintf(1, '%.4f; ', tp(i));
    end
end

fprintf(1, '\n\nDiscounted Pathogenic\n');
for i = size(tp, 2) : -1 : 1
    if i == 1
        fprintf(1, 'VS: ');
    elseif i == 2
        fprintf(1, 'ST: ');
    elseif i == 3
        fprintf(1, 'MO: ');
    else
        fprintf(1, 'SU: ');
    end
    
    if isnan(dtp(i))
        fprintf(1, '------; ');
    else
        fprintf(1, '%.4f; ', dtp(i));
    end
end

fprintf(1, '\n\n');

return



%[ThresholdB{10} flip(ThresholdP{10})]