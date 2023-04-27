function mask=mk_mask4Gspline(PT)
if(iscell(PT))
    cell_flg=1;
    P=PT{1};
    T=PT{2};
else
    cell_flg=0;
    P=PT(:,1);
    T=PT(:,2);
end
PTmask=[300000 14000
    269000 14000
    237000 14000
    232000 13000
    222000 12400
    213000 11200-1000
    204000 8700-700
    196000  7200-700
    180000  5400-700
    145000 4500-700
    97000  4300-700
    85000 4100-700
    73000 3000-300
    66000 2100-200
    62000 1300-200
    59000  873-100
    32000   720
    17000   720
    13000  700
    9800   630
    5400 440
    2500 320
    1400 280
    1000 245
    300 240
    273 235 %242  %235
    0 242];

if cell_flg
    nP=length(P);
    nT=length(T);
    mask=ones(nP,nT);
    [Pm,Tm]=ndgrid(P,T);
    for i=1:nP
        Tc=interp1(PTmask(:,1),PTmask(:,2),P(i));
        id= T<Tc;
        mask(i,id)=nan;
    end
    id= Pm<20 & Tm>500;
    mask(id)=nan;
    id= Pm<90 & Tm>500 & Tm<900;
    mask(id)=nan;
else
    nT=length(T);
    mask=ones(nT,1);
    for i=1:nT
        Tc=interp1(PTmask(:,1),PTmask(:,2),P(i));
        if(T(i)<Tc),mask(i)=nan;end
        if(P(i)<20 && T(i)>500),mask(i)=nan;end
        if(P(i)<90 && T(i)>500 && T(i)<900),mask(i)=nan;end
    end
end
    
    
    
    