function type = classifyInteraction(a1, a2, th)
    function c = getCat(a, t)
        if a > t
            c = 1;
        elseif a < -t
            c = -1;
        else
            c = 0;
        end
    end
    c1 = getCat(a1, th);
    c2 = getCat(a2, th);
    if c1 == 1 && c2 == 1
        type = 1;
    elseif (c1 == 1 && c2 == 0) || (c1 == 0 && c2 == 1)
        type = 2;
    elseif c1 == 0 && c2 == 0
        type = 3;
    elseif (c1 == -1 && c2 == 0) || (c1 == 0 && c2 == -1)
        type = 4;
    elseif (c1 == 1 && c2 == -1) || (c1 == -1 && c2 == 1)
        type = 5;
    elseif c1 == -1 && c2 == -1
        type = 6;
    else
        type = 3;
    end
end
