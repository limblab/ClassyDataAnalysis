function fitPPCA(binned)
%function fitPPCA
%       inputs:
%            binned data structure with pcaConfig filled with needed
%            parameters
%       Outputs:
%           puts data into ppcaData, and sets ranPPCAFit event
%           
%       binned.dimeReductionConfig Structure must have:
%
%           windows : start and end times of windows to perform
%           dimensionality reduction on (eg. [targetAppearsTrial1,
%           goCueTrial1;targetAppearsTrial2, goCueTrial2....]
%           
%           which : list of binned data column numbers to include in the 
%           dimensionality reduction
%
%           dimension: dimension of reduced dataset
%
%           segLength:

    
    [ppcaData.coeff,ppcaData.score,ppcaData.latent,ppcaData.mu,ppcaData.istropicVariance,ppcaData.stats]=ppca(binned.data{windows2mask(binned.data.t,binned.dimReductionConfig.windows),binned.dimReductionConfig.which},binned.dimReductionConfig.dimension);

    set(binned,'ppcaData', ppcaData);
    opData = binned.dimReductionConfig;
    evntData=loggingListenerEventData('fitPPCA',opData);
    notify(binned,'ranPPCAFit',evntData)
end

