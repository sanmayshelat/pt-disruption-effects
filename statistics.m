%% Statistics
totaltransfers = 0;
totalIVT = 0;
totalWT = 0;
for i = 1: numnodes %origin
    for j = 1: numnodes %destination
        for k = 1:numel(mastertransfers{i,j})
            totaltransfers = totaltransfers + ...
                mastertransfers{i,j}(k)*newload{i,j}(k);
            totalIVT = totalIVT + ...
                masterlength{i,j}(k)*newload{i,j}(k);
            totalWT = totalWT + ...
                masterwt{i,j}(k)*newload{i,j}(k);
        end
    end
end
totalutility = betaivt*totalIVT+betawt*totalWT+betatransfer*totaltransfers;
totaldemand = sum(sum(odmat));
averagetransfers = totaltransfers/totaldemand;
averageIVT = totalIVT/totaldemand;
averageWT = totalWT/totaldemand;
averageutility = totalutility/totaldemand;
            
            