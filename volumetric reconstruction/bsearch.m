function ind = bsearch(x, y)
    %function ind = bsearch(x, y)
    %find largest ind such that x(ind) <= y

    low = 1;
    high = length(x);
    while (low < high-1)
        ii = floor((high + low)/2);
        if (x(ii) > y)
            high = ii-1;
        else
            low = ii;
        end
    end
    if (x(high) <= y)
        ind = high;
    else
        ind = low;
    end
end