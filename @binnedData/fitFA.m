function fitFA(binned)
%function fitFA
%       inputs:
%            binned data structure with faConfig filled with needed
%            parameters
%       Outputs:
%           First n factors of data set specified in config file. Saved in
%           faData section of data structure
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

    dat = binned.dimRedHelper();
    kernSD = binned.dimReductionConfig.kernSD;
    runIdx =31;
    xDim = binned.dimReductionConfig.dimension;
    result = neuralTraj(runIdx,dat, 'method', 'fa', 'xDim', xDim, 'kernSDList', kernSD, 'segLength', binned.dimReductionConfig.segLength);
    faData = result;
    set(binned,'faData', faData);
    opData = binned.dimReductionConfig;
    evntData=loggingListenerEventData('fitFA',opData);
    notify(binned,'ranFAFit',evntData)
end

