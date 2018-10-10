for cellidx = 1:numel(def1)
    isSix = cellfun(@(x)isequal(x,def1{cellidx}),def);
    if isempty(find(isSix))
        disp(def1{cellidx})
    end
    
end