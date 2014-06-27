function finalOutput=code2func(xcode,inputList,outputList)
finalOutput=struct('code',cell2mat(xcode(:)'),'argins',...
    {inputList},'argouts',{outputList});
end
