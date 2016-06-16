function fitFA(binned)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 %this is a method function of the binnedData class and should be saved
    %in the @binnedData folder with the class definition and other methods
    %files
    %
    %bd.fitPds uses the configuration in the bd.pdConfig field to compute
    %preferred directions for each unit and stores the result in the
    %bd.pdData field
    %get our list of units
    method = 'fa';
    dat = dimRedHelper(binned, method);
    kernSD = 10;
    runIdx =20;
    xDim = binned.faConfig.dimension;
    result = neuralTraj(runIdx,dat, 'method', method, 'xDim', xDim, 'kernSDList', kernSD, 'segLength', binned.faConfig.segLength);
    faData = result;
    set(binned,'faData', faData);
    opData = binned.pcaConfig;
    evntData=loggingListenerEventData('fitFA',opData);
    notify(binned,'ranFAFit',evntData)
end

