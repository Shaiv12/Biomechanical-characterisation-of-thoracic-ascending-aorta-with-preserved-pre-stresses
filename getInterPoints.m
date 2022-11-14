%% getInterPoints

function d = getInterPoints(section1,section2)

numPoints = size(section1,1);

% d = zeros(2*numPoints);
for q = 1:1:numPoints  
    p1 = section1(q,:,:) ;
    p2 = section2(q,:,:);
    d.p{q} = [p1; p2];
end   

end