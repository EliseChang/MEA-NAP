Files = dir;
FileNames = cell(length(Files),1);
for k = 1:length(Files)
   name = erase(Files(k).name, ["'",".mat"]);
   FileNames{k,1} = name;
end
writecell(FileNames)