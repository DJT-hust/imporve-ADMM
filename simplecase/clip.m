function res=clip(x,xmin,xmax)

    if x<xmin
        res=xmin;
    elseif x>xmax
        res=xmax;
    else
        res=x;
    end
    