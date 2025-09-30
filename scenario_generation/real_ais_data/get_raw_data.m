function [data] = get_raw_data(source_dir,N)

fullPath = source_dir+"\traj_*.txt";
files = dir(char(fullPath));
n=length(files);
disp(files);

data = cell(N,1);

bar = waitbar(0,'Loading of data...');

if N>n 
    disp("Error: The number of tracks cannot exceed"+string(n));
else
    selectedInds = randsample(1:n,N);
    for i=1:N
        %disp(source_dir+"\test_"+set+string(i)+".csv");
        fileInd = selectedInds(i);
        fileName = string(files(fileInd).name);
        data{i} = readmatrix(source_dir+"\"+fileName);
        waitbar(i/N, bar);
    end
end

% if ~isempty(N)
%     if N>n 
%         disp("Error: The number of tracks cannot exceed"+string(n));
%     else
%         selectedInds = randsample(1:n,N);
%         for i=1:N
%             %disp(source_dir+"\test_"+set+string(i)+".csv");
%             fileInd = selectedInds(i);
%             fileName = string(files(fileInd).name);
%             data{i} = readmatrix(source_dir+"\"+fileName);
%             waitbar(i/N, bar);
%         end
%     end
% else
%     selectedInds = 1:n;
%         for i=1:n
%             %disp(source_dir+"\test_"+set+string(i)+".csv");
%             fileInd = selectedInds(i);
%             fileName = string(files(fileInd).name);
%             data{i} = readmatrix(source_dir+"\"+fileName);
%             waitbar(i/n, bar);
%         end

end

