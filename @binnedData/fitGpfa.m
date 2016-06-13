function fitGpfa(binned)
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
    method = 'gpfa';
    dat = dimRedHelper(binned, method);
    kernSD = 10;
    runIdx =100;
    xDim = binned.gpfaConfig.dimension;
    result = neuralTraj(runIdx,dat, 'method', method, 'xDim', xDim, 'kernSDList', kernSD, 'segLength', binned.gpfaConfig.segLength);
    gpfaData = result;
    set(binned,'gpfaData', gpfaData);
    opData = binned.gpfaConfig;
    evntData=loggingListenerEventData('fitGpfa',opData);
    notify(binned,'ranGPFAFit',evntData)
end

