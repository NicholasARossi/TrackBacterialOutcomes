%%% Snapshots ------------------------------------------

% movieprefix = 'marsal5-0';

% moviepos = '';

% moviepath = '../data/8May12/';

% 

% segpath = [moviepath, 'segmented_snaps', moviepos, '/'];

% 

% if ~exist(segpath, 'dir')

%   mkdir(segpath);

% end

 

%%% Movies ------------------------------------------

 

 

movieprefix = 'sokt';

moviepos = 'xy01';

moviepath = '/Users/nicholasrossi/Documents/work_docs/lab_work/2018/9_25_18/sokcfp/';

 
segpath = [moviepath, 'segmented_', moviepos, '/'];

 

if ~exist(segpath, 'dir')

  mkdir(segpath);

end

 

 

% Figure out the number of charachters associated with the frame number

% e.g.  movie_t01xy01c1.tif  is 2

%       movie_t01xy001c1.tif is 3

 

file_list=dir([moviepath, movieprefix, '*', moviepos, 'c1.tif']);

file1 = file_list(1).name;

 

frameN = length(file1) - length(movieprefix) - length(moviepos) - length('c1.tif');

 

