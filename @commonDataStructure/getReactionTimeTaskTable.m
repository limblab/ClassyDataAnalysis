function getReactionTimeTaskTable(cds,times)
    %this is a method function for the common_data_structure (cds) class, and
    %should be located in a folder '@common_data_structure' with the class
    %definition file and other method files
    %
    %computes the trial variables for the CO task and composes the trial
    %table in the cds using the task variables and the generic trial times
    %passed in from the calling function. This is intended to be called by 
    %the getTrialTable method of the cds class, rather than directly by a
    %user
    
    %get our word timing for changes in the state machine:
    % Isolate the individual word timestamps
    bumpWordBase = hex2dec('50');
    bumpMask=cds.words.word >= (bumpWordBase) & cds.words.word <= (bumpWordBase+5);
    bumpTimes = cds.words.ts(bumpMask)';
    bumpCodes = cds.words.word(bumpMask)';

    wordOTOn = hex2dec('40');
    otMask=bitand(hex2dec('f0'),cds.words.word) == wordOTOn;
    otOnTimes = cds.words.ts( otMask);
    otOnCodes = cds.words.word( otMask);
    
    wordGo = hex2dec('31');
    goCueTime = cds.words.ts(cds.words.word == wordGo);
    
    wordStim=hex2dec('60');
    stimMask=bitand(hex2dec('f0'),cds.words.word) == wordStim;
    stimTimes=cds.words.ts( stimMask );
    stimCodeList=cds.words.word( stimMask );
    
    %preallocate our trial variables:
    numTrials=numel(times.number);

    tgtSize=nan(numTrials,1);
    tgtRadius = nan(numTrials,1);
    tgtAngle=nan(numTrials,1);
    randomTargets = nan(numTrials,1);
    showTgtDuringBump = nan(numTrials,1);
    
    bumpTrial = nan(numTrials,1);
    bumpMagnitude=nan(numTrials,1);
    bumpAngle=nan(numTrials,1);
    bumpFloor = nan(numTrials,1);
    bumpCeiling = nan(numTrials,1);
    bumpStep = nan(numTrials,1);
    bumpRisePeriod=nan(numTrials,1);
    bumpHoldPeriod = nan(numTrials,1);
    
    stimTrial=false(numTrials,1);
    stimCode=nan(numTrials,1);
    
    isTrainingTrial=false(numTrials,1);
    
    recenterCursor=false(numTrials,1);
    hideCursor = nan(numTrials,1);
    
    intertrialPeriod=nan(numTrials,1);
    penaltyPeriod=nan(numTrials,1);
    ctrHold=nan(numTrials,1);
    bumpDelay=nan(numTrials,1);
    
    tgtOnTime = nan(numTrials,1);
    bumpTimeList = nan(numTrials,1);
    goCueList = nan(numTrials,1);
    
    %get the databurst version:
    dbVersion=cds.databursts.db(1,2);
    skipList=[];
    switch dbVersion
        case 1
            % loop thorugh our trials and build our list vectors:
            for trial = 1:numTrials
                %find and parse the current databurst:
                idxDB = find(cds.databursts.ts > times.startTime(trial) & cds.databursts.ts<times.endTime(trial), 1, 'first');
                if isempty(idxDB)
                    skipList=[skipList,trial];
                    continue
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
                %from mastercon code to ensure matching when extracting data from
                %databurst:
                % * Version 1 (0x01)
                % * ----------------
                % * byte  0:		uchar		=> number of bytes to be transmitted
                % * byte  1:		uchar		=> version number (in this case one)
                % * byte  2-3:	uchar		=> task code ('FC')
                % * bytes 4-5:	uchar       => version code
                % * byte  5-6:	uchar		=> version code (micro)
                % *
                % * bytes 7-10:  float		=> target angle
                % * byte  11:	uchar           => random target flag
                % * bytes 12-15: float		=> target radius
                % * bytes 16-19: float		=> target size
                % * byte  20:	uchar		=> show target during bump
                % *
                % * byte  21:                => bump trial flag
                % * bytes 22-25: float		=> bump direction
                % * bytes 26-29: float       => bump magnitude
                % * bytes 30-33: float		=> bump floor (minimum force(N) bump can take)
                % * bytes 34-37:	float		=> bump ceiling (maximum force(N) bump can take)
                % * bytes 38-41:	float		=> bump step
                % * bytes 42-45: float		=> bump duration
                % * bytes 46-49: float		=> bump ramp
                % *
                % * byte  50:	uchar		=> stim trial flag
                % * bytes 51:    uchar       => stim code
                % *
                % * byte  52:    uchar       => training trial flag
                % *
                % * byte  53:	uchar		=> recenter cursor flag
                % * byte  54:    uchar       => hide cursor during bump
                % *
                % * bytes 55-58: float		=> intertrial time
                % * bytes 59-62: float		=> penalty time
                % * bytes 63-67: float		=> bump hold time
                % * bytes 68-71: float		=> center hold time
                % * bytes 72-75: float		=> bump delay time
                % */
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                tgtAngle(trial)=bytes2float(cds.databursts.db(idxDB,7:10));
                randomTargets(trial)=cds.databursts.db(idxDB,11);
                tgtRadius(trial)=bytes2float(cds.databursts.db(12:15));
                tgtSize(trial)=bytes2float(cds.databursts.db(idxDB,16:19));
                showTgtDuringBump(trial) = cds.databursts.db(idxDB,20);

                bumpTrial(trial) = cds.databursts.db(idxDB,21);
                bumpAngle(trial)=bytes2float(cds.databursts.db(idxDB,22:25));
                bumpMagnitude(trial)=bytes2float(cds.databursts.db(idxDB,26:29));
                bumpFloor(trial) = bytes2float(cds.databursts.db(idxDB,30:33));
                bumpCeiling(trial) = bytes2float(cds.databursts.db(idxDB,34:37));
                bumpStep(trial) = bytes2float(cds.databursts.db(idxDB,38:41));
                bumpHoldPeriod(trial) = bytes2float(cds.databursts.db(idxDB,42:45));
                bumpRisePeriod(trial) = bytes2float(cds.databursts.db(idxDB,46:49));
                
                stimTrial(trial)= cds.databursts.db(idxDB,50);
                stimCode(trial) = cds.databursts.db(idxDB,51);
                
                isTrainingTrial(trial)=cds.databursts.db(idxDB,52);
                
                recenterCursor(trial)=cds.databursts.db(idxDB,53);
                hideCursor(trial)=cds.databursts.db(idxDB,54);

                intertrialPeriod(trial)=bytes2float(cds.databursts.db(idxDB,55:58));
                penaltyPeriod(trial)=bytes2float(cds.databursts.db(idxDB,59:62));
                ctrHold=bytes2float(cds.databursts.db(idxDB,67:70));
                bumpDelay=bytes2float(cds.databursts.db(idxDB,71:74));
                
                abortDuringBump = cds.databursts.db(idxDB,75);
                forceReaction = cds.databursts.db(idxDB,76);
                
                %now get things that rely only on words and word timing:
                idxOT=find(otOnTimes>times.startTime(trial) & otOnTimes < times.endTime(trial),1,'first');
                if isempty(idxOT)
                    tgtOnTime(trial)=nan;
                    %tgtID(trial)=nan; %target ID has no meaning in this version of the databurst
                else
                    tgtOnTime(trial)=otOnTimes(idxOT);
                    %tgtID(trial)=otOnCodes(idxOT); %target ID has no meaning in this version of the databurst
                end

                % Bump code and time
                idxBump = find(bumpTimes > times.startTime(trial) & bumpTimes < times.endTime(trial), 1, 'first');
                if isempty(idxBump)
                    bumpTimeList(trial) = nan;
                    %bumpList(trial) = nan;%bump ID has no meaning in this version of the databurst
                    bumpAngle(trial)=nan;
                else
                    bumpTimeList(trial) = bumpTimes(idxBump);
                    %bumpList(trial) = bitand(hex2dec('0f'),bumpCodes(idxBump));%bump ID has no meaning in this version of the databurst
                end

                % Go cue
                idxGo = find(goCueTime > times.startTime(trial) & goCueTime < times.endTime(trial), 1, 'first');
                if isempty(idxGo)
                    goCueList(trial) = nan;
                else
                    goCueList(trial) = goCueTime(idxGo);
                end

                %Stim code
                idx = find(stimTimes > times.startTime(trial) & stimTimes < times.endTime(trial),1,'first');
                if isempty(idx)
                    stimCode(trial) = nan;
                else
                    stimCode(trial) = bitand(hex2dec('0f'),stimCodeList(idx));%hex2dec('0f') is a bitwise mask for the trailing bit of the word
                end
            end

            %build table:
            trialsTable=table(ctrHold,tgtOnTime,goCueList,intertrialPeriod,penaltyPeriod,bumpDelay,bumpHoldPeriod,...
                                tgtSize,tgtAngle,tgtRadius,...
                                isTrainingTrial,...
                                bumpTimeList,abortDuringBump,bumpHoldPeriod,bumpRisePeriod,bumpMagnitude,bumpAngle,...
                                stimTrial,stimCode,... 
                                recenterCursor,forceReaction,hideCursor,...
                                'VariableNames',{'ctrHold','tgtOnTime','goCueTime','intertrialPeriod','penaltyPeriod',...
                                'bumpDelay','bumpHoldTime','tgtSize','tgtDir','tgtDistance',...
                                'isTrainingTrial',...
                                'bumpTime','abortDuringBump','bumpHoldPeriod','bumpRisePeriod','bumpMagnitude','bumpDir',...
                                'isStimTrial','stimCode',...
                                'recenterCursor','forceReaction','hideCursor'});

            trialsTable.Properties.VariableUnits={'s','s','s','s','s','s','s',...
                                                    'cm','deg','cm, cm','cm','bool','deg','deg','bool','pct',...
                                                    's','bool','bool','bool','bool','s','s','N','deg',...
                                                    'bool','pct','int',...
                                                    'bool','bool'};
            trialsTable.Properties.VariableDescriptions={'center hold time','outer target onset time','go cue time','intertrial time','penalty time','time after entering ctr tgt that bump happens','time after bump onset before go cue',...
                                                            'size of targets','angle of outer target','distance to outer target from center','only the correct target was shown',...
                                                            'time of bump onset','would we abort during bumps','the time the bump was held at peak amplitude',...
                                                            'the time the bump took to rise and fall from peak amplitude','magnitude of the bump','direction of the bump',...
                                                            'was there stimulation','code in the stim word',...
                                                            'did the cursor recenter after bump','did we force reaction time',,'did we hide the cursor'};
            
       
        otherwise
            error('getCObumpTaskTable:unrecognizedDBVersion',['the trial table code for BD is not implemented for databursts with version#:',num2str(dbVersion)])
    end
    
    trialsTable=[times,trialsTable];
    trialsTable.Properties.Description='Trial table for the CObump task';
    %sanitize trial table by masking off corrupt databursts with nan's:
    mask= ( trialsTable.ctrHold<0           | trialsTable.ctrHold>10000 | ...
            trialsTable.delayHold<0         | trialsTable.delayHold>10000 |...
            trialsTable.intertrialPeriod<0  | trialsTable.intertrialPeriod>10000 |...
            trialsTable.penaltyPeriod<0     | trialsTable.penaltyPeriod>10000 |...
            trialsTable.bumpHoldPeriod<0     | trialsTable.bumpHoldPeriod>10000 |...
            trialsTable.bumpRisePeriod<0     | trialsTable.bumpRisePeriod>10000 |...
            trialsTable.bumpMagnitude<-100     | trialsTable.bumpMagnitude>100 |...
            trialsTable.tgtSize<.000001);
    mask(skipList)=1;
    idx=find(mask);
    for j=5:size(trialsTable,2)
        if ~isempty(find(strcmp({'goCueTime','tgtOnTime','bumpTime','tgtID','bumpID'},trialsTable.Properties.VariableNames{j}),1))
            %skip things that are based on the words, not the databurst
            continue
        end
        if islogical(trialsTable{1,j})
            trialsTable{idx,j}=false;
        else
            trialsTable{idx,j}=nan(size(trialsTable{1,j}));
        end
    end
    
    set(cds,'trials',trialsTable)
    evntData=loggingListenerEventData('getCOTaskTable',[]);
    notify(cds,'ranOperation',evntData)
end