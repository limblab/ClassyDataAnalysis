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
    %The matching portions of this code are adapted from B-Dekleva's code
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
        

        unitsChannels=unique(units.data.chan);
        dataChannels=unique(data.chan);
        % [inBoth,unitsIdxList]=ismember(dataChannels,unitsChannels);
        disp(['compiling population statistics using method: ',method])
        %compute population level stats:
        if ~isempty(strfind(method,'shape') )
            %compute the population shape distribution:
            [ldaProj,coeff]=getShapeComps(units,data);
        end

        if ~isempty(strfind(units.appendConfig.method,'ISI') )
            %compute the population ISI distribution:
            [tsISI]=getISIcomps(units,data);
        end

        for i=1:numel(unitsChannels)
        %loop through all the channels in the existing units data
        disp(['working on channel: ',num2str(unitsChannels(i))])
            %find all units on this channel:
            unitsList=find(units.data.chan==unitsChannels(i));
            dataList=find(data.chan==unitsChannels(i));
            %clear/reset our working vars
            dataFlags=true(numel(dataList,1));
            unitFlags=true(numel(unitsList,1));
            tmpUnitData=[];
            %get invalid and unsorted and put them into tmpUnitData:
            [unitFlags,tmpUnitData]=addUnit(0,units.data,tmpUnitData,unitFlags,unitsList,0);
            [unitFlags,tmpUnitData]=addUnit(255,units.data,tmpUnitData,unitFlags,unitsList,0);
            [dataFlags,tmpUnitData]=addUnit(0,data,dataFlags,tmpUnitData,dataList,offset);
            [dataFlags,tmpUnitData]=addUnit(255,data,dataFlags,tmpUnitData,dataList,offset);
            %loop across all the units on this channel:
            for j=1:numel(unitsList)
                if ~unitFlags(j)
                    %if this unit in units was already handled, skip it
                    continue
                else
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
                                unitMean=mean(units(unitsList(j)).spikes.wave);
                                unitStdev=std(units(unitsList(j)).spikes.wave);
                                dataMean=mean(data(dataList(k).spikes.wave));
                                dataStdev=std(data(dataList(k).spikes.wave));

                                [alpha,alphaCI]=regress(unitMean,dataMean);
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

                            end
                        end
                    end
                end
                %use units.appendConfig.default to assign any unmatched units:
                tmpData=data(dataList(dataFlags));
                for k=1:numel(tmpData)
                    tmpData(k).spikes.ts=tmpData(k).spikes.ts+offset;
                end
                unmatched=[ units.data(unitsList(unitsFlags)) ; 
                           tmpData];
                       
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
                        %add our data to the unsorted unit
                        tmpUnitData(idx).spikes=sortrows([tmpUnitData(idx).spikes;unmatched(k)],'ts');
                    else
                        %append our data:
                        tmpUnitData=[tmpUnitData;unmatched(k)];
                    end
                end
                newUnitData=[newUnitData;tmpUnitData];
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
    uInfo.newUnitData.numUnits=numel(newUnitData);
    uInfo.newUnitData.numChan=numel(unique(newUnitData.chan));
    uInfo.newUnitData.hasSorting=~isempty(find([newUnitData.ID]>0 & [newUnitData.ID]<255,'first'));
    uInfo.sepTime=offset;
    uInfo.ksDist=tsISI;
    uInfo.shaprDist=ldaProj;
    evntData=loggingListenerEventData('appendData',uInfo);
    notify(units,'appended',evntData)
    
end
%local functions that don't share namespace with the main function:
function [flags,tmpUnitData]=addUnit(ID,UD,tmpUnitData,flags,uList,offset)
    chanIdx=find([UD(uList).ID]==ID);
    if ~isempty(chanIdx)
        tmpData=UD(uList(chanIdx)).spikes;
        tmpData.ts=tmpData.ts+offset;
        %check to see if we already have a matching unit:
        idx=find([tmpUnitData.ID]==ID);
        if ~isempty(idx)
            % merge if we do:
            tmpUnitData(idx).spikes=sortrows([tmpUnitData(idx).spikes;tmpData],'ts');
        else
            %just append the data
            tmpUnitData=[tmpUnitData;tmpData];
        end
        %now flag this unit as dealt with and move on 
        flags(chanIdx)=false;
    end
end
function [ldaProj,coeff]=getShapeComps(units,data)
    %loop through all units and get mean and stdev of wave
    numUnits=numel(units.data);
    unitsMean=nan(numUnits,size(units.data(1).spikes.wave,2));
    unitsStdev=unitsMean;
    unitsCount=nan(numUnits,1);
    for i=1:numUnits
        unitsMean(i,:)=mean(units.data(i).spikes.wave);
        unitsStdev(i,:)=std(units.data(i).spikes.wave);
        unitsCount(i)=size(units.data(i).spikes,1);
    end
    numData=numel(data);
    dataMean=nan(numData,size(data(1).spikes.wave,2));
    dataStdev=dataMean;
    dataCount=nan(numUnits,1);
    for i=1:numData
        dataMean(i,:)=mean(data(i).spikes.wave);
        dataStdev(i,:)=std(data(i).spikes.wave);
        dataCount(i)=size(data(i).spikes,1);
    end
    %concatenate data together:
    allMean=[unitsMean;dataMean];
    allStdev=[unitsStdev;dataStdev];
    allCount=[unitsCount;dataCount];
    chans=[units.data.chan;data.chan];
% % % %     %now loop through the units and construct d' distances for all units on
% % % %     %different channels:
% % % %     sepMat=nan(0.5*numel(chans)^2,size(allMean,2)+1);
% % % %     sameChan=nan(size(sepMat,1));
% % % %     idx=1;
% % % %     for i=1:numel(chans)
% % % %         for j=i:numel(chans)
% % % %             sameChan= chans(i)==chans(j);
% % % %             %get scaling factor:
% % % %             [alpha,alphaCI]=regress(allMean(i,:),allMean(j,:));
% % % %             %convert CI into stdev:
% % % %             alphaStdev=diff(alphaCI)/(2*1.96);
% % % %             %now get distance between the waves in wavespace:
% % % %             dist=alpha*allMean(i,:)-allMean(j,:);
% % % %             jointStdev=sqrt(allStdev(i)+allStdev(j));
% % % %             %now use stdev and values to put d' measures into sepMat:
% % % %             sepMat(idx,:)=[dist./jointStdev , alpha/alphaStdev];
% % % %             idx=idx+1;
% % % %         end
% % % %     end
% % % %     %now clear any the residual nans so they don't cause problems later.
% % % %     %Also use abs to convert to separation, rather than signed distance:
% % % %     sepMat=abs(sepMat(isnan(sepMat(:,1)),:));
    numWaves=size(allMean,1);
    numPoints=size(allMean,2);
    waveMat=repmat(reshape(allMean,[1,numWaves,numPoints]),[numWaves,1,1]);
    stdevMat=repmat(reshape(allStdev,[1,numWaves,numPoints]),[numWaves,1,1]);
    spikesMat=repmat(allCount,[1,numWaves,numPoints]);
    alphaMat=nan(numWaves,numWaves);
    alphaStdMat=alphaMat;
    for i=1:size(waveMat,1)
        for j=1:size(waveMat,2)
            %now get the alpha (gain) factor to multiply the transpose
            %element in order to match the scale of the non-transpose
            %element
            [alphaMat(i,j),alphaCI]=regress(squeeze(waveMat(i,j,:)),squeeze(waveMat(j,i,:)));
            alphaStdMat(i,j)=diff(alphaCI)/(2*1.96);
        end
    end
    %compute all distances in wavespace. Subtract the scaled transpose from
    %the non-transpose (recall that alphas are scales to apply to the
    %transpose). note this results in a symmetric matrix, with duplicate 
    %entries.
    waveMat=waveMat-permute(waveMat,[2,1,3])*repmat(alphaMat,[1,1,numPoints]);
    %create an upper triangular mask excluding the self comparisons and 
    %duplicates in the bottom half
    mask=triu(true(numWaves),1);
    %use the mask to get a list of differences in wavespace
    diffs=waveMat(repmat(mask,[1,1,numPoints]));
    %reshape the differences into a column matrix where each row is a
    %difference observation include the alpha values here as the last 
    %difference:
    diffs=[reshape(diffs,numel(diffs)/numPoints,numPoints), alphaMat(mask)];
    %now compute joint standard deviation (S1*N1+S2*N2)/(N1+N2):
    stdevMat=sqrt((stdevMat.*spikesMat+permute(stdevMat,[2,1,3]).*permute(spikesMat,[2,1,3]))  ./ ...
                    (spikesMat+permute(spikesMat,[2,1,3])));
    %convert the 3D standard deviation matrix into a 2D matrix to match the
    %diffs matrix. Again, we add on the value for the alpha
    stdevs=stdevMat(repmat(mask,[1,1,numPoints]));
    stdevs=[reshape(stdevs,numel(stdevs)/numPoints,numPoints),alphaStdMat(mask)];
    %calculted dPrime from the differences and standard deviations
    dPrime=diffs./stdevs;
    %now get the logical index for things that are the same channel:
    diffMask=true(numWaves,numWaves);
    for i=1:numel(chans)
        tmp=find(allChans==chans(i));
        diffMask(tmp,tmp)=false;
    end
    knownDiff=diffMask(mask);
    
    %now lets make an LDA classifier to maximally separate units on same
    %and different channels:
    [~,~,~,~,coeff]=classify([],dPrime,knownDiff);
    %now project our same channel data onto the LDA axis:
    ldaProj = dPrime(knownDiff,:)*coeff(1,2).linear;
    ldaProj=sort(ldaProj);
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