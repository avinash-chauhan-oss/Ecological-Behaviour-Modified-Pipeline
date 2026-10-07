function alpha = calcAlpha(gc, gm)
    if isnan(gm) || gm < 1e-6
        if isnan(gc) || gc < 1e-6
            alpha = 0.0;
        else
            alpha = 1.0;
        end
    else
        alpha = (gc - gm) / gm;
    end
end
