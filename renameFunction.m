function parameterName = renameFunction(parameterName)

parameterName(parameterName == '!') = '';
parameterName(parameterName == '-') = '';
parameterName(parameterName == '/') = '_';
if parameterName(length(parameterName)) == '.'
    parameterName(length(parameterName)) = '';
end

if parameterName(length(parameterName)) == ')'
    parameterName(length(parameterName)) = '';
end


specialCharecterIndex = strfind(parameterName,' ');
for var1 = 1:1:length(specialCharecterIndex)
    newParameterName = strrep(parameterName(specialCharecterIndex(var1) + 1), ... 
        parameterName(specialCharecterIndex(var1) + 1),upper(parameterName(specialCharecterIndex(var1) + 1)));
    parameterName(specialCharecterIndex(var1) + 1) = newParameterName;
end

specialCharecterIndex = strfind(parameterName,'(');
for var1 = 1:1:length(specialCharecterIndex)
    newParameterName = strrep(parameterName(specialCharecterIndex(var1) + 1), ... 
        parameterName(specialCharecterIndex(var1) + 1),upper(parameterName(specialCharecterIndex(var1) + 1)));
    parameterName(specialCharecterIndex(var1) + 1) = newParameterName;
end

specialCharecterIndex = strfind(parameterName,')');
for var1 = 1:1:length(specialCharecterIndex)
    newParameterName = strrep(parameterName(specialCharecterIndex(var1) + 1), ... 
        parameterName(specialCharecterIndex(var1) + 1),upper(parameterName(specialCharecterIndex(var1) + 1)));
    parameterName(specialCharecterIndex(var1) + 1) = newParameterName;
end

specialCharecterIndex = strfind(parameterName,'.');
for var1 = 1:1:length(specialCharecterIndex)
    newParameterName = strrep(parameterName(specialCharecterIndex(var1) + 1), ... 
        parameterName(specialCharecterIndex(var1) + 1),upper(parameterName(specialCharecterIndex(var1) + 1)));
    parameterName(specialCharecterIndex(var1) + 1) = newParameterName;
end

parameterName(parameterName == ' ') = '';
parameterName(parameterName == '(') = '';
parameterName(parameterName == ')') = '';
parameterName(parameterName == '.') = '';


end