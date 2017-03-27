function [testLLforMaxValid, maxValid, maxValidInd, maxParams, maxTestLL_TODO_REMOVE, maxParams_testTODO_REMOVE] = findBestParams(name )

testLLforMaxValid=[]; maxValid=[]; maxValidInd=[]; maxParams=[];

% beginningString = 'RunAll_results_';

nameLen = length(name);

%%find all files starting with beginningString
resFiles = dir(['*']);

%%load all results from files
% trainLL = nan(1,nRes);
% validLL = nan(1,nRes);
% testLL = nan(1,nRes);

trainLL = [];
validLL = [];
testLL = [];

res = {};
params = {};
for i=1:length(resFiles)
    %read params
    file = fopen(resFiles(i).name,'r');
    if file==-1
        continue;
    end
    
    fgetl(file);
    
    %find row starting with name
    dataLine=-1;
    ind = 0;
    while true
        ind = ind+1;
        line = fgetl(file);
        if isnumeric(line) && line==-1
            break;
        end
        if length(line)>=nameLen && strcmp(line(1:nameLen),name)
            dataLine = ind;
            break;
        end
    end
    fclose(file);
    if (dataLine==-1) %%relevant line not found
        continue;
    end
    
    %read data
    lastRes = [];
    try
        lastRes = dlmread(resFiles(i).name,'\t',[dataLine,1,dataLine,3]);
    catch exception
        fprintf('WARNING: could not read data in %s\n',resFiles(i).name);
        continue;
    end
    %     trainLL(:,i) = res{i}(:,1);
    %     validLL(:,i) = res{i}(:,2);
    %     testLL(:,i) = res{i}(:,3);
    
    assert(~isempty(lastRes));
    
    %append results and params
    res{end+1} = lastRes;
    params{end+1} = fileread(resFiles(i).name);
    trainLL = [trainLL, lastRes(:,1)];
    validLL = [validLL, lastRes(:,2)];
    testLL = [testLL, lastRes(:,3)];
end

if (~isempty(res))
    %%find max for validation LL, for all experiments
    [maxValid, maxValidInd] = max(validLL,[],2);
    [maxTestLL_TODO_REMOVE, maxTest_ind] = max(testLL,[],2); %should not be used as test performance (otherwise, overfitting) but just as check
    [testLLforMaxValid] = testLL(maxValidInd);
    
    %find best parameters
    maxParams = params{maxValidInd};
    maxParams_testTODO_REMOVE = params{maxTest_ind};
else
    fprintf('No scores found for dataset %s \n',name);
    testLLforMaxValid = nan;
end


end