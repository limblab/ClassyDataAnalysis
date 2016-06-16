function fitGpfa(binned)
%function fitGPFA
%       inputs:
%            binned data structure with gpfaConfig filled with needed
%            parameters
%       Outputs:
%           First n factors of data set specified in config file. Saved in
%           gpfaData section of data structure
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
    method = 'gpfa';
    dat = dimRedHelper(binned, method);
    kernSD = binned.gpfaConfig.kernSD;
    runIdx =100;
    xDim = binned.gpfaConfig.dimension;
    result = neuralTraj(runIdx,dat, 'method', method, 'xDim', xDim, 'kernSDList', kernSD, 'segLength', binned.gpfaConfig.segLength);
    gpfaData = result;
    set(binned,'gpfaData', gpfaData);
    opData = binned.gpfaConfig;
    evntData=loggingListenerEventData('fitGpfa',opData);
    notify(binned,'ranGPFAFit',evntData)
end

