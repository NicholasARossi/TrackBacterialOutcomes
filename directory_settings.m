


movieprefix = 'ppit';

moviepos = 'xy10';

moviepath = '/Users/nicholasrossi/Documents/work_docs/lab_work/2018/9_25_18/ppi/';

 
segpath = [moviepath, 'segmented_', moviepos, '/'];

 

if ~exist(segpath, 'dir')

  mkdir(segpath);

end

 

file_list=dir([moviepath, movieprefix, '*', moviepos, 'c1.tif']);

file1 = file_list(1).name;

 

frameN = length(file1) - length(movieprefix) - length(moviepos) - length('c1.tif');
 
