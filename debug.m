clear all, close all
movieprefix = 'sokt';

moviepos = 'xy01';

moviepath = '/Users/nicholasrossi/Documents/work_docs/lab_work/2018/9_25_18/sokcfp/';

 
segpath = [moviepath, 'segmented_', moviepos, '/'];

 


posmod=moviepos;
posmod(regexp(posmod,'[0]'))=[];
%%% ingredients of clist data, data3D, def def3D gate neighbor
load([segpath,posmod,'/clist.mat']);




% 
load([segpath,posmod,'/seg/', movieprefix, '01', moviepos,'_seg','.mat'],'mask_cell');

pos=data(check_number,[18,19]);
Imout=bwlabel(mask_cell);
validate1=Imout;
validate1([validate1~=check_number])=0;
% 
% RGB = insertMarker(validate1,pos);
imshow(validate1);
