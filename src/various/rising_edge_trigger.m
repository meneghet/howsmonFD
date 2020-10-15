function ynew = rising_edge_trigger(y)
% Apply rising edge trigger

ynew = zeros(length(y),1);

x = find(y);
if ~isempty(x)
    edges_= [x(1); x(find(diff(x) > 1) + 1)];
    ynew(edges_) = 1;
end

end