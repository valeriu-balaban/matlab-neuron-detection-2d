function Points = merge_random(Points, merge_dis)

times = ones(size(Points,1),1);
temp = 0;
while true
    num = size(Points, 1);
    if temp == num
        break;
    end
    temp = num;
    for i = 1:num-1
        for j = i+1:num
            dis = (Points(i, :) - Points(j, :)).^2;
            dis = sqrt(sum(dis));
            if dis < merge_dis
                Points(j,:) = NaN;
                times(j) = NaN;
            end
        end
    end
    Points(isnan(Points(:,1)), :) = [];
end