function fitPPCA(binned)
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
    method = 'ppca';
    dat = dimRedHelper(binned, method);
    kernSD = 10;
    runIdx =20;
    xDim = binned.ppcaConfig.dimension;
    result = neuralTraj(runIdx,dat, 'method', method, 'xDim', xDim, 'kernSDList', kernSD, 'segLength', binned.gpfaConfig.segLength);
    ppcaData = result;
    set(binned,'ppcaData', ppcaData);
    opData = binned.ppcaConfig;
    evntData=loggingListenerEventData('fitPPCA',opData);
    notify(binned,'ranPPCAFit',evntData)
end

