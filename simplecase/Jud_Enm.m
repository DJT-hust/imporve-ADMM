function jud=Jud_Enm(OD,E_nm,epsilon)

    jud=0;
    N=size(OD,1);
    for n=1:N
        if abs(E_nm(OD(n,1),OD(n,2))+E_nm(OD(n,2),OD(n,1)))>epsilon
            jud=1;
            break;
        end
    end