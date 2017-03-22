function fitPPCA(binned)
%function fitPPCA
%       inputs:
%            binned data structure with ppcaConfig filled with needed
%            parameters
%       Outputs:
%           First n factors of data set specified in config file. Saved in
%           ppcaData section of data structure
%       Config Structure:
%
%           windows : start and end times of windows to perform
%           dimensionality reduction on (eg. [targetAppearsTrial1,
%           goCueTrial1;targetAppearsTrial2, goCueTrial2....]
%           
%           dimension: dimension of reduced dataset
%
%           segLength: Lenght of smaller segments to cut trials into. Makes
%           computation of factors faster for FA and GPFA if trials are of
%           equal length
%
%           trials: trial numbers corresponding to windows
%
%           kernSD: smoothing kernel standard deviation. Larger value acts
%           as a larger low pass filter on trajectories
    method = 'ppca';
    dat = dimRedHelper(binned, method);
    kernSD = binned.ppcaConfig.kernSD;
    runIdx =20;
    xDim = binned.ppcaConfig.dimension;
    result = neuralTraj(runIdx,dat, 'method', method, 'xDim', xDim, 'kernSDList', kernSD, 'segLength', binned.gpfaConfig.segLength);
    ppcaData = result;
    set(binned,'ppcaData', ppcaData);
    opData = binned.ppcaConfig;
    evntData=loggingListenerEventData('fitPPCA',opData);
    notify(binned,'ranPPCAFit',evntData)
end

