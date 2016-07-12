function getUNTTaskTable(cds,times)
    %this is a method of the commonDataStructure class and must be saved in
    %the @commonDataStructure folder with the other method definitions
    %
    %produces a trials table with the following information:
    %startTime      -time of trial start
    %ctrOnTime      -time the center target appeared
    %tgtOnTime      -time the ourter target appeared
    %goCueTime      -time of the go cue
    %endTime        -time the trial ended
    %result         -result of the trial (RAFI)
    %tgtDir         -direction of the true target in degrees
    %cuePrior       -kappe for the prior of the cue
    %tgtPrior       -kappa for the prior of the actual targets
    

    centerOnTime       = cds.words.ts(cds.words.word==hex2dec('30'));
    otOnTime        = cds.words.ts(cds.words.word==hex2dec('40'));
    goCueTime       = cds.words.ts(cds.words.word==hex2dec('31'));
    %preallocate vectors:
    numTrials=numel(times.number);
    tgtDirList=nan(numTrials,1);
    tgtPriorList=nan(numTrials,1);
    cuePriorList=nan(numTrials,1);
    goTime=nan(numTrials,1);
    OTTime=nan(numTrials,1);
    ctrOnTime=nan(numTrials,1);
    % For each trial complete code
    for trial=1:numTrials
        %find the databurst for this trial
        idxDB = find(cds.databursts.ts > times.startTime(trial) & cds.databursts.ts<times.endTime(trial), 1, 'first');
        %get target and prior info from databurst
        if ~isempty(idxDB) 
            tgtDirList(trial) = 180*bytes2float(cds.databursts.db(idxDB,10:13))/pi;
            tgtPriorList(trial) = bytes2float(cds.databursts.db(idxDB,14:17));
            cuePriorList(trial)=bytes2float(cds.databursts.db(idxDB,18:21));
            if cuePriorList > 100000 
                cuePriorList(trial) = NaN;
            end
        else
            tgtDirList(trial) = NaN;
            tgtPriorList(trial) = NaN;
            cuePriorList(trial) = NaN;
        end
        % get the timestamp for the go cue
        gT = goCueTime(find(goCueTime<times.endTime(trial) & goCueTime>times.startTime(trial),1,'first'));
        if isempty(gT)
            goTime(trial)=NaN;
        else
            goTime(trial)=gT;
        end
        % get the timestamp for the outer target appearance 
        OTT = otOnTime(find(otOnTime<times.endTime(trial) & otOnTime>times.startTime(trial),1,'first'));
        if isempty(OTT)
            OTTime(trial)=NaN;
        else
            OTTime(trial)=OTT;
        end
        %get the timestamp for center target appearance:
        cOT = find(centerOnTime<times.endTime(trial) & centerOnTime>times.startTime(trial),1,'first');
        if isempty(cOT)
            ctrOnTime(trial)=NaN;
        else
            ctrOnTime(trial)=cOT;
        end
    end

    % Deal with weird prior databursts
    checkPrior = @(burst) burst < 10e-5 | burst > 1e5+1 | isnan(burst);
    badBursts = find(checkPrior(tgtPriorList));
    goodBursts = find(~checkPrior(tgtPriorList));
    for i = 1:length(badBursts)
        bb = badBursts(i);
        ind_dists = abs(goodBursts - bb);
        replacer_ind = goodBursts(find(ind_dists==min(ind_dists),1,'first'));
        replacer = tgtPriorList(replacer_ind);
        tgtPriorList(bb)= replacer;
    end

    trialsTable=table(roundTime(ctrOnTime,.001),roundTime(goTime,.001),roundTime(OTTime,.001),...
                        tgtDirList,cuePriorList,tgtPriorList,...
                        'VariableNames',{'ctrOnTime','goCueTime','tgtOnTime','tgtDir','cuePrior','tgtPrior'});
    trialsTable.Properties.VariableUnits={'s','s','s','Deg','AU','AU'};
    trialsTable.Properties.VariableDescriptions={'center target onset time','go cue time','outer target onset time',...
                                                    'actual target direction','kappa of the von mises function for the visual cue','kappa for the von mises function for the actual target locations'};
    
    trialsTable=[times,trialsTable];
    trialsTable.Properties.Description='Trial table for the UNT task';
    set(cds,'trials',trialsTable)
    evntData=loggingListenerEventData('getCOTaskTable',[]);
    notify(cds,'ranOperation',evntData)
end