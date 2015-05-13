function fastq_dump_in_parallel( sra_dir , processors)
% FASTQ_DUMP_IN_PARALLEL fastq-dump is slow, so uncompress files up in parallel
% SRA_dir: directory of SRA files
% processors: number of processors to run in parallel

matlabpool(processors);

cd(sra_dir);
files = dir('.');

unix('chmod 666 *.sra')

file_names = {files.name};
rem_inds = [1 1 zeros(1,numel(file_names)-2)];
for i=3:numel(file_names)
   if isempty(findstr('.sra', file_names{i}))
      rem_inds(i) = 1;
   end
end
file_names(rem_inds==1) = [];

for i=1:numel(file_names)
   base_names{i} = file_names{i}(1:end-4);
end

parfor i=1:numel(file_names)
   unix(['fastq-dump ', sra_dir, '/', file_names{i}])
   unix(['rm ', file_names{i}])
   unix(['gzip ', base_names{i}, '.fastq'])
end

matlabpool close

end

