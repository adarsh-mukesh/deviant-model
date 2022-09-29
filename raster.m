function[]=raster(fire_time)
for i=1:size(fire_time,1)
    arr=fire_time(i,:);
    arr=arr(arr~=0);
    for j=1:numel(arr)
        line([arr(j) arr(j)],[i-0.4 i+0.4],'color','r');
    end
end
xlabel('Time(ms)');
ylabel('neuron number');
end
        
    