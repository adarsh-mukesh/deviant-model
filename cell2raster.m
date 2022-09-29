function []=cell2raster(fire_time,clr)
for i=1:length(fire_time)
    arr=fire_time{i,:};
    if isempty(arr)
        continue
    end
    arr=fix(arr*1000);
    for j=1:numel(arr)
        scatter([arr(j) arr(j)],[i i],4,'MarkerEdgeColor',clr,'MarkerFaceColor',clr)
        hold on
    end
end
end