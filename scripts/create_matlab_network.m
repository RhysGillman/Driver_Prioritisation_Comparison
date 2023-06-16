function [] = main_SCS( network_file )

% read in the csv file
data = readtable(network_file, 'ReadVariableNames', true, 'ReadRowNames', false);

%remove the unwanted column
data = data(:, {'protein_1','protein_2'});

% extract the unique gene names from both protein_1 and protein_2 columns
node0 = unique([data.protein_1; data.protein_2]);

% create the edge0 cell array
edge0 = cell(height(data), 2);
edge0(:,1) = data.protein_1;
edge0(:,2) = data.protein_2;

% save the variables to a .mat file
save("../tmp/tmp_network.mat", "edge0", "node0");

end