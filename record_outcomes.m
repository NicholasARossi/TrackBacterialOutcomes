function record_outcomes(fr_start,fr_end)

%% Known issues/bugs/code updated needed
% * only works for one color (channel c2)
% * need to click in window before you can record an outcome
% * would be better if you could crop/zoom
% * little red outline in subsequent frames sort of sucks without expanding
% image

% I needed to initialize the following to avoid static workspace errors.
% The ones that are needed get assigned later.
movieprefix = []; moviepos = []; moviepath = []; segpath = []; Imout = []; file_list = []; file1 = []; frameN = [];

directory_settings; % get path info
close all

disp('-------------------------------------------------------------------')
disp(' Record Outcomes Commands:')
disp(' ')
disp(' a = mark cell as lysed')
disp(' s = mark cell as non-growing')
disp(' d = mark cell as filamented')
disp(' f = mark cell as exited the frame')
disp(' ')
disp(' u = undo last assignment')
disp(' ')
disp('-------------------------------------------------------------------')
disp(' ')

outcomes = {'a', 's', 'd', 'f'};
outcomes_names = {'lysed', 'non-growing', 'filamented', 'exited'};


fr_end = length(dir([moviepath, movieprefix, '*', moviepos, 'c1.tif'])); % find the number of phase image files

% Read in all images and store in Imall
for fr = fr_start:fr_end
    clear Imorig rim pim Im
    Imorig = imread([moviepath, movieprefix, strN(fr, frameN), moviepos, 'c1.tif']);
    
    % This loads the PI channel
    rim = imread([moviepath, movieprefix,  strN(fr, frameN), moviepos, 'c3.tif']);

    Im = double(Imorig);
    Im = Im/maxmax(Im);
    
    rscale = [4000 16000];
    
    rim(rim>rscale(2))=rscale(2);
    rim(rim<rscale(1))=rscale(1);

    rim = (rim - rscale(1))/(rscale(2)-rscale(1));
    
    Imallorig{fr} = Im+double(rim);
end


% If the cells in fr_start have not already been identified with
% locate_cells, then find them. Otherwise use what already exists.
if ~exist([segpath,moviepos,'/seg/', movieprefix, '01', moviepos,'_seg','.mat']) > 0
    for n =1:3
        copyfile([moviepath, movieprefix, '01', moviepos, 'c',num2str(n),'.tif'],[segpath, movieprefix, '01', moviepos, 'c',num2str(n),'.tif'])
    end
    %runs the supper segger analysis of the first frame
    processExp(segpath)
end
list = {'marRAB','acrAB','InaA','SoxS','ompF','hdeA','purA','micF','gadX','sodA','tolC','crp','sigma70' ...
'lacuv5','rpsT','Fis','H-NS','rrnbp1','rob','lacZ','tonB','CodB','ompC','dnaQ'};
[indx,tf] = listdlg('PromptString','Select Promoter:','ListString',list);

Experiment={list{indx},'50 carb'};
tempvar=strsplit(moviepos,'y');
posmod=moviepos;

if str2num(tempvar{2})<10
    posmod(regexp(posmod,'[0]'))=[];
end
    
load([segpath,posmod,'/seg/', movieprefix, '01', moviepos,'_seg','.mat'],'mask_cell') % Load segmentation data for fr_start
load([segpath,posmod,'/clist.mat'],'data', 'data3D', 'def', 'def3D');

[arrayRows,arrayCols] = size(data);
NewCol = NaN(arrayRows,1);
%Add new column
data = [data NewCol];
def={def{:},'Time of Death'};

Imout=bwlabel(mask_cell);

% Set up GUI with slider
hFig = figure('menu','none');
hAx = axes('Parent',hFig);
hSlider = uicontrol('Parent',hFig, 'Style','slider', 'Value', fr_start, 'Min', fr_start,...
    'Max', fr_end, 'SliderStep', [1 5]/(fr_end-fr_start), ...
    'Position', [150 5 300 20], 'Callback', @slider_callback);
hTxt = uicontrol('Style','text', 'Position',[290 28 20 15], 'String', num2str(fr_start));


% Callback function
    function slider_callback(hObj, ~)
        fr = round(get(hObj,'Value'));      % get frame from slider
        imshow(Imall{fr}, 'Parent',hAx)     % show new image
        set(hTxt, 'String',num2str(fr))     % update text on slider
    end




cell = 1;
last_outcome = '';

while cell <= max(max(Imout)) % all labeled cells need to be assigned an outcome
    
    disp(['Cell #: ', num2str(cell)])
    
    % Color code cells on first image
    for i = 1:3
        Im3(:,:,i) = Imallorig{fr_start};
    end
    [idx, idy] = find(Imout < cell & Imout > 0); for i = 1:length(idx); Im3(idx(i), idy(i), 1:3) = [0.5 0.5 0.7]; end % all mark cells colored gray-blue
    [idx, idy] = find(Imout == cell);            for i = 1:length(idx); Im3(idx(i), idy(i), 1:3) = [1 0 0]; end % current cell colored red
    Imall{fr_start} = Im3;
    imshow(Imall{fr_start}, 'Parent', hAx);
    set(hSlider, 'Value', fr_start)
    set(hTxt, 'String', num2str(fr_start))
    
    % Mark current cell position in all subsequent frames
    mask = zeros(size(Imout));
    mask(Imout == cell) = 1;
    bwp = bwperim(mask,8);
    [idx, idy] = find(bwp);
    for fr = fr_start+1:fr_end
        for i = 1:3
            Im3(:,:,i) = Imallorig{fr};
        end
        for i = 1:length(idx)
            Im3(idx(i), idy(i), 1:3) = [1 0 0];
        end
        Imall{fr} = Im3;
    end
    
    disp('Record the outcome of the cell indicated in red.')
    w = waitforbuttonpress;
    while w == 0
        w = waitforbuttonpress;
    end
    cc=get(gcf,'currentcharacter');
    disp(cc)
    if isequal(cc, 'u') % undo last assignment, assuming their was one
        cell = cell - 1;

    elseif sum(strcmp(cc, outcomes)) > 0 % first make sure that the selected charachter is one of the possible outcomes  
        for oi = 1:length(outcomes)
            if isequal(cc, outcomes{oi})
                if isequal(cc, 'a')
                    data(cell,end)=fr;
                elseif isequal(cc,'d')
                    data(cell,end)=-1;
                elseif isequal(cc,'s')
                    data(cell,end)=0;
                elseif isequal(cc,'f')
                    data(cell,end)=NaN(1,1);
                end
            end
        end

        cell = cell + 1;
    else
        disp(['"', cc, '" is not a valid selection. Please try again.'])
        last_outcome = '';
    end
    
end



% Load the fr_start image data for each fluor




save([segpath,posmod,'/clist.mat'], 'data', 'data3D', 'def', 'def3D','Experiment');
[y,Fs] = audioread('finish.wav');
sound(y,Fs);
end