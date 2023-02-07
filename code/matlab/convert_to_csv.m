% Ecoli
expInd = 2;
tbl = array2table(data{expInd,11}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'ecoli_exometabolites.csv')
tbl = array2table(data{expInd,12}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'ecoli_exometabolites_std.csv')

% OD
tbl = table(data{expInd,2}, data{expInd,5}, data{expInd,6}, 'VariableNames', {'Time', 'OD mean', 'OD std'});
writetable(tbl, 'ecoli_OD.csv')

% Glucose
tbl = table(data{expInd,2}, data{expInd,3}, data{expInd,4}, 'VariableNames', {'Time', 'Glucose mean', 'Glucose std'});
writetable(tbl, 'ecoli_glucose.csv')

%%
% B licheniformis WT batch
expInd = 3;
tbl = array2table(data{expInd,11}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'b_licheniformis_exometabolites.csv')
tbl = array2table(data{expInd,12}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'b_licheniformis_exometabolites_std.csv')

% OD
tbl = table(data{expInd,2}, data{expInd,5}, data{expInd,6}, 'VariableNames', {'Time', 'OD mean', 'OD std'});
writetable(tbl, 'b_licheniformis_OD.csv')

% Glucose
tbl = table(data{expInd,2}, data{expInd,3}, data{expInd,4}, 'VariableNames', {'Time', 'Glucose mean', 'Glucose std'});
writetable(tbl, 'b_licheniformis_glucose.csv')

%%
% Yeast
expInd = 4;
tbl = array2table(data{expInd,11}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'yeast_exometabolites.csv')
tbl = array2table(data{expInd,12}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'yeast_exometabolites_std.csv')

% OD
tbl = table(data{expInd,2}, data{expInd,5}, data{expInd,6}, 'VariableNames', {'Time', 'OD mean', 'OD std'});
writetable(tbl, 'yeast_OD.csv')

% Glucose
tbl = table(data{expInd,2}, data{expInd,3}, data{expInd,4}, 'VariableNames', {'Time', 'Glucose mean', 'Glucose std'});
writetable(tbl, 'yeast_glucose.csv')

%%
% C glutamicum
expInd = 10;
tbl = array2table(data{expInd,11}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'c_glutamicum_exometabolites.csv')
tbl = array2table(data{expInd,12}, 'VariableNames', data{expInd, 10});
tbl2 = [ table(data{expInd,2}, 'VariableNames', {'Time'})  tbl];
writetable(tbl2, 'c_glutamicum_exometabolites_std.csv')

% OD
tbl = table(data{expInd,2}, data{expInd,5}, data{expInd,6}, 'VariableNames', {'Time', 'OD mean', 'OD std'});
writetable(tbl, 'c_glutamicum_OD.csv')

% Glucose
tbl = table(data{expInd,2}, data{expInd,3}, data{expInd,4}, 'VariableNames', {'Time', 'Glucose mean', 'Glucose std'});
writetable(tbl, 'c_glutamicum_glucose.csv')