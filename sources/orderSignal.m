function [se,ord_ind] = orderSignal(s,se)

[num_sources,~] = size(s);
ord_ind = zeros(1,num_sources);

for i=1:num_sources
    max = 1e-16;
    for j=1:num_sources
        
        %         val = sum(mean(abs(s(i,:)-se(j,:)).^2));
        val = abs(sum(diag(flipud(corrcoef(s(i,:),se(j,:)))))/2);
        if val > max;
            max = val;
            ind_max = j;
        end
    end
    ord_ind(1,i) = ind_max;
end
se=se(ord_ind,:);

end

