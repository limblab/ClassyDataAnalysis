function appendData(units,data,varargin)
    %this is a method function of the unitData class and should be found in
    %the @unitData folder.
    %
    %accepts the units field from a cds and appends the data onto
    %the current units structure. appendUnits uses the value in
    %offset to time shift all the spikes of the new data
    %
    %appendData accesses the appendConfig property of the unitData class
    %in order to append new units onto the existing unit data. This
    %property defines the tests and thresholds that appendData will use to
    %determine whether units are the same across multipel files. The
    %appendConfig proptery MUST have the following fields:
    %method:    'shape' | 'ISI' | 'shapeISI' | 'number'
    %           -shape uses cohen's D (d') between pairs of waves, and an
    %           LDA classifier to estimate whether a pair of waves is drawn
    %           from the distribution of 2 different units, or the
    %           distribution of same units
    %           -ISI uses the KS-test statistic to parameterize differences
    %           in ISI between 2 units. the ISI parameter should only be
    %           used for merging the same task, as task differences WILL
    %           change the ISI properties of a unit.
    %           -shapeISI uses both the shape and ISI as described above.
    %           P-values for the 2 tests are multiplied together and
    %           compared to threshold^2
    %           -number just matches the unit number and assumes that the
    %           user has sorted the files so that the units are the same.
    %threshold: # on range 0-1 corresponding to p-value for rejecting
    %           similiarity
    %default:   'unsorted' | 'invalid' | 'delete'
    %           -unsorted relabels unmatched units as unsorted
    %           -invalid relabels unmatched units as invalid
    %           -delete removes unmatched units from unitData
    %
    %appendData does not return anything, instead, the unitData.data field
    %is updated. appendData also returns information via a notification to
    %listeners, including the number of units in the new data, as well as
    %the distributions of shape and ISI parameters in the dataset. The
    %experiment class has a listener that will catch this notification and
    %store the data, in the operationLog property (inherited from the
    %operationLogger class)
    %
    %The unit matching portions of this code are adapted from B-Dekleva's code
    %that is based on the method outlined in Rebesco etal, 'Rewiring neural
    %interactions by micro-stimulation' in frontiers in systems
    %neuroscience from 2012. The major modification here is that we use
    %Cohen's D (d') rather than simple euclidean norms to compare wave
    %shapes, thus capturing the noise within a specific wave as part of
    %our discrimination. The Rabesco result relies on the variance of
    %means, ignoring the variation of waves within a single unit
    %

    if nargin>3
        error('appendData:tooManyInputs','appendData accpts up to 3 imputs')
    end

    if ~isempty(varargin)
        offset=varargin{1};
    else
        offset=[];
    end
    newUnitData=[];
    tsISI=[];
    ldaProj=[];
    %if units is empty, simply fill it
    if isempty(units.data)
        if ~isempty(offset) && offset>0
            %this case is here to make the unit times behave the
            %same way as timeSeriesData in the case where the
            %user passes an offset. In theory this case should
            %never be used
            warning('appendData:shiftedNewData','applying a time shift to data that is being placed in an empty unitData.data field')
            for i=1:length(data)
                data(i).spikes.t=data(i).spikes.t+offset;
            end
        end

        set(units,'data',data)
    else
        %sanity checks:
        %do we have the same arrays?
        diffArrays=setdiff({units.data.array},{data.array});
        if ~isempty(diffArrays)
            error('appendUnits:differentArrays',['this unitData has the following array(s): ',strjoin(unique({units.data.array}),','), ' while the new units structure has the following array(s): ',strjoin(unique({data.array}),',')])
        end

        %do we have the same unit set?
        diffUnits=setdiff([cell2mat({units.data.chan}),cell2mat({units.data.ID})],[cell2mat({data.chan}),cell2mat({data.ID})]);
        if ~isempty(diffUnits)
            warning('appendUnits:differentUnits',['the new units field has ',num2str(numel(diffUnits)),'different units from the units in this unitData structure'])                
        end
        %is offset larger than the biggest value in the original
        %data?
        f=@(x) max(x.ts);
        maxUnitsTime=max(cellfun(f,{units.data.spikes}));
        if isempty(offset)
            offset=maxUnitsTime+1;
        end
        if maxUnitsTime>offset
            error('appendUnits:inadequateOffset','The offset for timestamps must be larger than the maximum timestamp in the existing data. Suggest using the duration of existing timeseries data like kinematics to estimate an offset.');
        end
        %ok we passed the sanity checks, now update the time of all
        %the spikes by adding offset, and append data to units

                                %         %build tables with the channel, unitID, and array to
                                %         %server as unique keys for the old and new unitdata
                                %         unitsKey=table([units.data.chan]',[units.data.ID]',char({units.data.array}'),'VariableNames',{'chan','ID','array'});
                                %         dataKey=table([data.chan]',[data.ID]',char({data.array}'),'VariableNames',{'chan','ID','array'});
                                %         %now handle the stuff in dataKey that's in unitsKey:
                                %         [inBoth,inBothIdx]=ismember(dataKey,unitsKey);
                                %         %directly assign elements of units.data, rather than
                                %             %using set so we can avoid copying the whole units.data
                                %             %field and wasting a bunch of memory. This is still
                                %             %really slow, but I can't figure out how to correct it
                                %             %given the structure of our units data.
                                %         for i=1:length(inBoth)
                                %             if inBoth(i)
                                %                 data(i).spikes.ts=data(i).spikes.ts+offset;
                                %                 units.data(inBothIdx(i)).spikes=[ units.data(inBothIdx(i)).spikes ; data(i).spikes ];
                                %             end
                                %         end
                                %         %now handle the stuff that's only in the dataKey
                                %         inDataOnly=find(~inBoth);
                                %         if ~isempty(inDataOnly)
                                %             units.data(end+1:end+length(inDataOnly))=data(inDataOnly);
                                %         end
        method=units.appendConfig.method;
        newUnitData=[];
        

        unitsChannels=unique([units.data.chan]);
        arrays=unique({units.data.array});
        dataChannels=unique([data.chan]);
        % [inBoth,unitsIdxList]=ismember(dataChannels,unitsChannels);
        disp(['compiling population statistics using method: ',method])
        %compute population level stats:
        if ~isempty(strfind(method,'shape') )
            %compute the population shape distribution:
            [ldaProj,coeff]=getShapeComps(units,data,units.appendConfig.SNRThreshold);
        end

        if ~isempty(strfind(units.appendConfig.method,'ISI') )
            %compute the population ISI distribution:
            [tsISI]=getISIcomps(units,data);
        end
        %handle each array separately since we can duplicate channel
        %numbers across arrays:
        for m=1:numel(arrays)
            currArr=arrays{m};
            for i=1:numel(unitsChannels)
                %loop through all the channels in the existing units data
                disp(['working on array: ',currArr,' channel: ',num2str(unitsChannels(i))])
                %find all units on this array&channel:
                unitsList=find([units.data.chan]==unitsChannels(i) & strcmp({units.data.array},currArr));
                dataList=find([data.chan]==unitsChannels(i) & strcmp({data.array},currArr));
                %clear/reset our working vars
                dataFlags=true(numel(dataList),1);
                unitFlags=true(numel(unitsList),1);
                tmpUnitData=[];
                %get invalid and unsorted and put them into tmpUnitData:
                [unitFlags,tmpUnitData]=addUnit(0,units.data,tmpUnitData,unitFlags,unitsList,0);
                [unitFlags,tmpUnitData]=addUnit(255,units.data,tmpUnitData,unitFlags,unitsList,0);
                [dataFlags,tmpUnitData]=addUnit(0,data,tmpUnitData,dataFlags,dataList,offset);
                [dataFlags,tmpUnitData]=addUnit(255,data,tmpUnitData,dataFlags,dataList,offset);
                
                
                
                %loop across all the units on this channel:
                for j=1:numel(unitsList)
                    if ~unitFlags(j)
                        %if this unit in units was already handled, skip it
                        continue
                    else
                    unitMean=mean(units.data(unitsList(j)).spikes.wave);
                    unitStdev=std(units.data(unitsList(j)).spikes.wave);
                    %skip this unit if the SNR is too low
                    %compute SNR as range of mean wave normalized by mean
                    %stdev of each point.
                    SNR=(max(unitMean)-min(unitMean))/mean(unitStdev);
                    if SNR<units.appendConfig.SNRThreshold;
                        continue
                    end
                    
                    %test all units in data for match, and append
                    %them if match is found

                        for k=1:numel(dataList)
                            if ~dataFlags(k) || ~strcmp(units.data(unitsList(j)).array,data(dataList(k)).array)
                                %if this unit in data was already handled, or
                                %if the array name doesn't match, skip it
                                continue
                            else
                                %compare data(dataList(k)) to
                                %units.data(unitsList(j)) :

                                if ~isempty(strfind( method,'shape'))
                                    %get the difference in the mean waveshapes\
                                    dataMean=mean(data(dataList(k)).spikes.wave);
                                    dataStdev=std(data(dataList(k)).spikes.wave);
                                    %skip this unit if the SNR is too low
                                    %compute SNR as range of mean wave normalized by mean
                                    %stdev of each point.
                                    SNR=(max(dataMean)-min(dataMean))/mean(dataStdev);
                                    if SNR<units.appendConfig.SNRThreshold;
                                        continue
                                    end
                                    [alpha,alphaCI]=regress(unitMean',dataMean');
                                    alphaStdev=diff(alphaCI)/(2*1.96);

                                    dPrime=[abs(unitMean-dataMean)./sqrt(dataStdev+unitStdev) , alpha/alphaStdev];
                                    %now project dPrime onto the LDA axis:
                                    proj = dPrime*coeff(1,2).linear;
                                    %now see where we are on the distribution of
                                    %known different cells to estimate p-value:
                                    pShape=sum(ldaProj>=proj)/numel(ldaProj);
                                else
                                    pShape=1;
                                end
                                if ~isempty(strfind(method, 'ISI'))
                                    [~,~,ks]=kstest2(units.data(unitsList(j)).spikes.ts,data(dataList(k)).spikes.ts);
                                    %now see where we are on the distribution of
                                    %known different cells to estimate p-value:
                                    pISI=min([sum(tsISI>=ks),sum(tsISI<=ks)]/numel(tsISI))/2;%the min and divide by 2 business is to handle the 2sided nature of the statistic
                                else
                                    pISI=1;
                                end
                                if ~isempty(strfind(method, 'number'))
                                    %check to see if we have a number match, and if set
                                    %an arbitrary pval that is larger than the
                                    %threshold so we trigger a merge, otherwise set one
                                    %smaller than the threshold so we skip merging
                                    if units.data(unitsList(j)).ID==data(dataList(k)).ID
                                        pNum=units.appendConfig.threshold+1;
                                    else
                                        pNum=units.appendConfig.threshold-1;
                                    end
                                else
                                    pNum=1;
                                end
                                %now merge p-vals:
                                %at this point we need to think about how we will compare
                                %to threshold in cases where we have multiple p-values. for
                                %instance if we have threshold=.05, and we have 2 tests, do
                                %we want both to attain significance? if one is SUPER
                                %significant is that enough? In this instance we are
                                %asserting that one highly significant value can compensate
                                %for a marginal/non-significant p-value. We will insist
                                %that the joint p-value for 2 merged tests be equivalent to
                                %threshold^2. That is if threshold=.05, then the joint p
                                %must be less than .05*.05=.0025 to attain significance.
                                %Since we don't want to mess with things later, we will
                                %simply compute the joint P here and then sqrt it to let us
                                %compare to the original threshold value:
                                pVal=(pShape*pISI*pNum)^(1/sum([~isempty(strfind( method,'shape')),~isempty(strfind( method,'ISI')),~isempty(strfind( method,'number'))]));
                                if pVal>=units.appendConfig.threshold
                                    %merge
                                    tmpUnitData=[tmpUnitData;units.data(unitsList(j))];
                                    tmpData=data(dataList(k)).spikes;
                                    tmpData.ts=tmpData.ts+offset;
                                    tmpUnitData(end).spikes=[tmpUnitData(end).spikes;tmpData];
                                    unitFlags(j)=false;
                                    dataFlags(k)=false;
                                    %since we found a match, break from the
                                    %loop across the units in data
                                    break
                                end
                            end
                        end
                    end
                    %use units.appendConfig.default to assign any unmatched units:
                    tmpData=data(dataList(dataFlags));
                    for k=1:numel(tmpData)
                        tmpData(k).spikes.ts=tmpData(k).spikes.ts+offset;
                    end
                    unmatched=[ units.data(unitsList(unitFlags)) ; 
                               tmpData'];

                    for k=1:numel(unmatched)
                        switch units.appendConfig.default
                            case 'unsorted'
                                idx=find([tmpUnitData.ID]==0);
                            case 'invalid'
                                idx=find([tmpUnitData.ID]==255);
                            case 'delete'
                                %just skip adding unmatched units to newUnitData
                                idx=[];
                            otherwise
                                error('appendData:badDefault',['appendData does not recognize the option: ', units.appendConfig.default,'. please update unitData.set:appendConfig so that it catches this error, or update this method to handle this option'])
                        end
                        if ~isempty(idx)
                            %add our data to the unsorted/invalid unit
                            tmpUnitData(idx).spikes=sortrows([tmpUnitData(idx).spikes;unmatched(k).spikes],'ts');
                        else
                            %append our data:
                            tmpUnitData=[tmpUnitData;unmatched(k)];
                        end
                    end
                    newUnitData=[newUnitData;tmpUnitData];
                end

            end
        end
        %now that we have handled every channel, put our data into units:
        set(units,'data',newUnitData)
    end
    uInfo.added.numUnits=numel([data.ID]);
    uInfo.added.numChan=numel(unique([data.chan]));
    uInfo.added.hasSorting=~isempty(find([data.ID]>0 & [data.ID]<255,1,'first'));
    uInfo.inUnits.numUnits=numel([units.data.ID]);
    uInfo.inUnits.numChan=numel(unique([units.data.chan]));
    uInfo.inUnits.hasSorting=~isempty(find([units.data.ID]>0 & [units.data.ID]<255,1,'first'));
    uInfo.sepTime=offset;
    uInfo.distribution.ksDist=tsISI;
    uInfo.distribution.shaprDist=ldaProj;
    evntData=loggingListenerEventData('appendData',uInfo);
    notify(units,'appended',evntData)
    
end
%local functions that don't share namespace with the main function:
function [flags,tmpUnitData]=addUnit(ID,UD,tmpUnitData,flags,uList,offset)
    chanIdx=find([UD(uList).ID]==ID);
    if ~isempty(chanIdx)
        tmpData=UD(uList(chanIdx));
        tmpData.spikes.ts=tmpData.spikes.ts+offset;
        %check to see if we already have a matching unit:
        idx=[];
        if ~isempty(tmpUnitData)
            idx=find([tmpUnitData.ID]==ID);
        end
        if ~isempty(idx)
            % merge if we do:
            tmpUnitData(idx).spikes=sortrows([tmpUnitData(idx).spikes;tmpData.spikes],'ts');
        else
            %just append the whole struct to tmpUnitData
            tmpUnitData=[tmpUnitData;tmpData];
        end
        %now flag this unit as dealt with and move on 
        flags(chanIdx)=false;
    end
end
function [ldaProj,coeff]=getShapeComps(units,data,SNRThresh)
    %loop through all units and get mean and stdev of wave
    disp('compiling sorted units already in unitData')
    numUnits=numel(units.data);
    unitsMean=nan(numUnits,size(units.data(1).spikes.wave,2));
    unitsStdev=unitsMean;
    unitsCount=nan(numUnits,1);
    unitsChans=[units.data.chan];
    unitsArray={units.data.array};
    unitsInUnits=ones(size(unitsChans));
    for i=1:numUnits
        if(units.data(i).ID==0 || units.data(i).ID==255)
            %don't bother with invalid or unsorted
            continue
        end
        unitsMean(i,:)=mean(units.data(i).spikes.wave);
        unitsStdev(i,:)=std(units.data(i).spikes.wave);
        unitsCount(i)=size(units.data(i).spikes,1);
    end
%     range=max(unitsMean,[],2)-min(unitsMean,[],2);
%     SNR=range./mean(unitsStdev,2);
%     %build mask to remove undesireable units:
%     mask=(~isnan(unitsCount) & ... stuff leftover from units 0 and 255
%             unitsChans<128 & ... sorted units on the analog front panel of the cerebus
%             SNR>=SNRThresh);% lowSNR units
    
    %clear out the empty rows from unsorted and invalid units
    mask=~isnan(unitsCount);
    unitsMean=unitsMean(mask,:);
    unitsStdev=unitsStdev(mask,:);
    unitsCount=unitsCount(mask);
    unitsChans=unitsChans(mask);
    unitsArray=unitsArray(mask);
    unitsInUnits=unitsInUnits(mask);
    %check for and clear out any sorted stuff from the analog front panel:
    mask=unitsChans<128;
    if ~isempty(find(~mask))
        warning('appendData:unitsHasHighChannel','the existing unit data has sorted units above 128 channels. This is not normally neural data so we are removing it. If you need this data, you should take care of it prior to merging units')
        unitsMean=unitsMean(mask,:);
        unitsStdev=unitsStdev(mask,:);
        unitsCount=unitsCount(mask);
        unitsChans=unitsChans(mask);
        unitsArray=unitsArray(mask);
        unitsInUnits=unitsInUnits(mask);
    end
    %now remove anything with max SNR below SNRThresh
    range=max(unitsMean,[],2)-min(unitsMean,[],2);
    SNR=range./mean(unitsStdev,2);
    mask=SNR>=SNRThresh;
    unitsMean=unitsMean(mask,:);
    unitsStdev=unitsStdev(mask,:);
    unitsCount=unitsCount(mask);
    unitsChans=unitsChans(mask);
    unitsArray=unitsArray(mask);
    unitsInUnits=unitsInUnits(mask);
    %now work on the data
    disp('compiling sorted units already in new data')
    numData=numel(data);
    dataMean=nan(numData,size(data(1).spikes.wave,2));
    dataStdev=dataMean;
    dataCount=nan(numUnits,1);
    dataChans=[data.chan];
    dataArray={data.array};
    dataInUnits=zeros(size(dataChans));
    for i=1:numData
        if(data(i).ID==0 || data(i).ID==255)
            continue
        end
        dataMean(i,:)=mean(data(i).spikes.wave);
        dataStdev(i,:)=std(data(i).spikes.wave);
        dataCount(i)=size(data(i).spikes,1);
    end
    %clear out the empty rows from unsorted and invalid units
    mask=~isnan(dataCount);
    dataMean=dataMean(mask,:);
    dataStdev=dataStdev(mask,:);
    dataCount=dataCount(mask);
    dataChans=dataChans(mask);
    dataArray=dataArray(mask);
    dataInUnits=dataInUnits(mask);
    %check for and clear out any sorted stuff from the analog front panel:
    mask=dataChans<128;
    if ~isempty(find(~mask))
        warning('appendData:unitsHasHighChannel','the existing unit data has sorted units above 128 channels. This is not normally neural data so we are removing it. If you need this data, you should take care of it prior to merging units')
        dataMean=dataMean(mask,:);
        dataStdev=dataStdev(mask,:);
        dataCount=dataCount(mask);
        dataChans=dataChans(mask);
        dataArray=dataArray(mask);
        dataInUnits=dataInUnits(mask);
    end
    %now remove anything with max SNR below SNRThresh
    range=max(dataMean,[],2)-min(dataMean,[],2);
    SNR=range./mean(dataStdev,2);
    mask=SNR>=SNRThresh;
    dataMean=dataMean(mask,:);
    dataStdev=dataStdev(mask,:);
    dataCount=dataCount(mask);
    dataChans=dataChans(mask);
    dataArray=dataArray(mask);
    dataInUnits=dataInUnits(mask);
    %concatenate data& units together:
    disp('merging units')    
    allMean=[unitsMean;dataMean];
    allStdev=[unitsStdev;dataStdev];
    allCount=[unitsCount;dataCount];
    allChans=[unitsChans';dataChans'];
    allArray=[unitsArray';dataArray'];
    allInUnits=[unitsInUnits';dataInUnits'];
    
    %now loop through the units and get scaling factors:
    numWaves=size(allMean,1);
    numPoints=size(allMean,2);
    waveMat=repmat(reshape(allMean,[1,numWaves,numPoints]),[numWaves,1,1]);
    stdevMat=repmat(reshape(allStdev,[1,numWaves,numPoints]),[numWaves,1,1]);
    spikesMat=repmat(allCount,[1,numWaves,numPoints]);
    disp('computing scaling factors for best shape-matching')
    alphaMat=nan(numWaves,numWaves);
    alphaCIMat=nan(numWaves,numWaves,2);
    for i=1:size(waveMat,1)
        for j=1:size(waveMat,2)
            %now get the alpha (gain) factor to multiply the transpose
            %element in order to match the scale of the non-transpose
            %element
            [alphaMat(i,j),alphaCIMat(i,j,:)]=regress(squeeze(waveMat(i,j,:)),squeeze(waveMat(j,i,:)));
%            alphaStdMat(i,j)=diff(alphaCI)/(2*1.96);
        end
    end
    %compute all distances in wavespace. Subtract the scaled transpose from
    %the non-transpose (recall that alphas are scales to apply to the
    %transpose). note this results in a symmetric matrix, with duplicate 
    %entries.
    disp('computing distances metrics')
    waveMat=waveMat-permute(waveMat,[2,1,3]).*repmat(alphaMat,[1,1,numPoints]);
    %create an upper triangular mask excluding the self comparisons and 
    %duplicates in the bottom half
    mask=triu(true(numWaves),1);
    %use the mask to get a list of differences in wavespace
    diffs=abs(waveMat(repmat(mask,[1,1,numPoints])));
    %reshape the differences into a column matrix where each row is a
    %difference observation include the alpha values here as the last 
    %difference:
    %convert alpha into 1 sided distribution:
    alphaMat=abs(alphaMat-1);
    diffs=[reshape(diffs,numel(diffs)/numPoints,numPoints), alphaMat(mask)];
    %diffs=[sqrt(sum(reshape(diffs,numel(diffs)/numPoints,numPoints).^2,2)), alphas];
    %now compute joint standard deviation (S1*N1+S2*N2)/(N1+N2):
    stdevMat=stdevMat.*spikesMat;
    stdevMat=stdevMat+permute(stdevMat,[2,1,3]);
    stdevMat= stdevMat./ (spikesMat+permute(spikesMat,[2,1,3]));
    %convert the 3D standard deviation matrix into a 2D matrix to match the
    %diffs matrix. Again, we add on the value for the alpha
    stdevs=stdevMat(repmat(mask,[1,1,numPoints]));
    alphaStdMat=diff(alphaCIMat,1,3);
    %convert CI into stdev for the scaling factor. Remember to convert CI
    %values for alphas<1 so the CI range matches the range for the
    %converted alpha:
    alphaStdMat=abs(alphaStdMat)/(2*1.96);
    stdevs=[reshape(stdevs,numel(stdevs)/numPoints,numPoints),alphaStdMat(mask)];
    %stdevs=[mean(reshape(stdevs,numel(stdevs)/numPoints,numPoints),2),alphaStdMat(mask)];
    %calculte dPrime from the differences and standard deviations. Log
    %transform to get from a positive only, skewed distribution to 
    %something that looks normal
    dPrime=log(diffs./stdevs);
    %now get the logical index for things that are the same channel:
    disp('getting projections onto LDA')
    diffMask=true(numWaves,numWaves);
    for i=1:numel(allChans)
        tmp=find(allChans==allChans(i) & strcmp(allArray,allArray{i}) & allInUnits~=allInUnits(i));
        diffMask(i,tmp)=false;
        diffMask(tmp,i)=false;
    end
    knownDiff=diffMask(mask);
    %now lets get the axis between the mean cluster position for our dPrime
    %data. This is the same as the LDA axis, but we get to skip all the
    %logic associated with classifying individual points:
    coeff=mean(dPrime(knownDiff,:))-mean(dPrime(~knownDiff,:));
    ldaProj=sort(dPrime(knownDiff,:)*coeff');
end
function [tsISI]=getISIcomps(units,data)
    numUnits=numel(units.data);
    numData=numel(data);
    tsISI=nan(0.5*(numUnits+numData)^2,1);
    idx=1;
    for i=1:numel(chans)
        for j=i:numel(chans)
            if i<=numUnits
                if j<=numUnits
                    [~,~,tsISI(idx)]=kstest2(units(i).spikes.ts,units(j).spikes.ts);
                else
                    [~,~,tsISI(idx)]=kstest2(units(i).spikes.ts,data(j-numUnits).spikes.ts);
                end
            else
                if j<=numUnits
                    [~,~,tsISI(idx)]=kstest2(data(i-numUnits).spikes.ts,units(j).spikes.ts);
                else
                    [~,~,tsISI(idx)]=kstest2(data(i-numUnits).spikes.ts,data(j-numUnits).spikes.ts);
                end
            end
            idx=idx+1;
        end
    end
    tsISI=sort(tsISI);
    numISI=numel(tsISI);
end