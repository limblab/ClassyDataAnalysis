function [ dat ] = dimRedHelper( binned, method )
    switch method
        case 'gpfa'
            windows = binned.gpfaConfig.windows;
            for i = 1:length(windows)
                dat(i).trialId = binned.gpfaConfig.trials(i);
                dat(i).spikes(:,:) = table2array(binned.bins(find(binned.bins.t>windows(i,1),1):find(binned.bins.t>windows(i,2),1), 12:end))'./1000;
            end
        case 'pca'
            windows = binned.pcaConfig.windows;
            for i = 1:length(windows)
                dat(i).trialId = binned.pcaConfig.trials(i);
                dat(i).spikes(:,:) = table2array(binned.data(find(binned.data.t>windows(i,1),1):find(binned.data.t>windows(i,2),1), 12:end))'./1000;
            end
        case 'ppca'
            windows = binned.ppcaConfig.windows;
            for i = 1:length(windows)
                dat(i).trialId = binned.ppcaConfig.trials(i);
                dat(i).spikes(:,:) = table2array(binned.bins(find(binned.bins.t>windows(i,1),1):find(binned.bins.t>windows(i,2),1), 12:end))'./1000;
            end
        case 'fa'
            windows = binned.faConfig.windows;
            for i = 1:length(windows)
                dat(i).trialId = binned.faConfig.trials(i);
                dat(i).spikes(:,:) = table2array(binned.bins(find(binned.bins.t>windows(i,1),1):find(binned.bins.t>windows(i,2),1), 12:end))'./1000;
            end
        otherwise
            error('not an approved method. Shame on you')
    end
end

