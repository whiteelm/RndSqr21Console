    program Test5
    implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
    dimension ztemp(58), betam(58), qwork(2000)
    dimension z1(58), z2(58), z3(58), z(58), x2(58), x3(58), x4(58)
    dimension dC(4),dL(4), Um(4), aC(9), bC(4), dZ0(4), er(2), dCC(4), dLL(4), em(2)
    a = 16;
    b = 13.;
    d1 = 3.; d2 = 3.;
    h1 = 0.65; h2 = 0.65;
    s=2.9; s01 = 3.55;
    e = 1;
    s01 = s01 + d1/2
    s = s + d2/2 + d1/2
    s02 = a - s01 - s
    pi =  dacos(-1.d0)		  
    data zi   /(0.d0, 1.d0) / 
    n = 2      
    nn = 58    
    r1 = d1/2.
    r2 = d2/2.
    iprint = -2
    iguess = 1
    m=24;
    z1(1) = s01
    do i = 1, m+1
        alf = 2./m
        z1(i+1) = zi*(h1+r1) + r1*exp(-zi*(i+5)*alf*pi) + s01;
    end do
    z1(27) = s01
    z1(28) = s01 + s
    do i = 1, m+1
        z1(i+28)= zi*(h2+r2) + r2*exp(-zi*(i+5)*alf*pi) + s01 + s;
    end do
    z1(54) = s01 + s
    z1(55) = s01 + s + s02
    z1(56) = s01 + s + s02 + b*zi
    z1(57) = b*zi
    z1(58) = 0
    betam(1) = -0.5
    betam(2) = -0.5 + alf/2.
    do i = 1, m-1
        betam(i+2) = alf;
    end do
    betam(26) = -0.5 + alf/2.
    betam(27) = -0.5
    betam(28) = -0.5
    betam(29) = -0.5 + alf/2.
    do i = 1, m-1
        betam(i+29) = alf;
    end do
    betam(53) = -0.5 + alf/2.
    betam(54) = -0.5
    betam(55) = -0.5
    betam(56) = -0.5
    betam(57) = -0.5
    betam(58) = -0.5
    z1c = dcmplx(b/2, a/2)
    nptsq = 6
    call qinit(nn,betam,nptsq,qwork)
    do 1 k = 1,nn
        z(k) = exp(dcmplx(0.d0, k-nn))
1   continue
    tol = 1d-6
    call scsolv(iprint,iguess,tol,errest,nn,c,z,z1c,z1,betam,nptsq,qwork)
    z20 =zi;				   
    z201=-zi;					
    do 2 k = 1, nn
        z2(k)=(z(k)*z201-z20)/(z(k)-1)  
        x2(k)=dreal(z2(k));         
2   continue
    do 3 k = 1,nn-4
        x3(k) = 2*(x2(k) - x2(2)) / (x2(nn-5) - x2(2)) - 1.
3   continue
    x4(1)=x3(1);
    x4(2)=x3(2);
    x4(3)=x3(26);
    x4(4)=x3(27);
    x4(5)=x3(28);
    x4(6)=x3(29);
    x4(7)=x3(53);
    x4(8)=x3(54);
    n2 = 3
    M=1000
    call GHIONE(x4,aC,n2,M)
    call refor(aC,bC,n)
    print*,'Capacitance matrix all in air bC'
    call dprint(bC,n)
    call capa(bC,dC,e,n)
    dC = dC * e
    print*,'Capacitance matrix [C] (pF/m)'
    call dprint(dC,n)
    call induc(bC,dL,n)
    print*,'Inductance matrix [L] (nH/m)'
    call dprint(dL,n)
    dCC=dC 
    dLL=dL
    call dminv(dLL,n,ad)
    call nroot(n,dCC,11.127*dLL,em,Um)
    print*,'Modal voltages matrix [Um] (V)'
    call dprint(Um,n)
    print*,'Er [em]'
    do 33 k = 1,n
33  print*,em(k)
    call impedance(n,dC,Um,em,dZ0)
    print*,'Impedance matrix [Z0] (Ohm)'
    call dprint(dZ0,n)
    end

