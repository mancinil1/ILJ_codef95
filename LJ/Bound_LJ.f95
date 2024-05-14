program BSwave_1D
implicit none
real*8 :: r0_L,r0_R,Psi0_L,Psi0_R,Psi1_L,Psi1_R,E_tot,error
real*8 :: tol,d1_L,D1_R,sumR,E1,n_integL,n_integR,sumL
real*8 :: norm_coef,areaL,areaR,area,rescale,dE,Ei,En,E2
real*8, allocatable :: M_L(:,:),M_R(:,:)
integer :: i,n_tot,n,qn
print *,'Starting initial E calculation.'
print *,'==============================================================================='
 tol=0.00001d0		!Check
 dE=1.0d0/10.0d0**5	!Check
 open(11,file='input1.txt')
 read(11,*)n_tot,E_tot,Psi0_L,Psi0_R,r0_L,r0_R,qn 
 close(11)
 allocate(M_L(n_tot+2,2),M_R(n_tot+2,2))
 do i=0,n_tot
    Ei=E_tot+(i*dE)
    Psi1_R=1.0d0/10.0d0**10	!Check
    Psi1_L=1.0d0/10.0d0**10	!Check
    call numerov(n_tot,Ei,Psi0_L,Psi1_L,r0_L,M_L)
    r0_L=M_L(1,1)
    Psi0_L=M_L(1,2)
    call numerov(n_tot,Ei,Psi0_R,Psi1_R,r0_R,M_R)
    r0_R=M_R(1,1)
    Psi0_R=M_R(1,2)
    rescale=M_L(n_tot+1,2)/M_R(n_tot+1,2)
    do n=1,n_tot+2
	M_R(n,2)=M_R(n,2)*rescale
    end do      
    d1_L=(M_L(n_tot+2,2)-M_L(n_tot+1,2))/(M_L(n_tot+2,1)-M_L(n_tot+1,1))
    d1_R=(M_R(n_tot,2)-M_R(n_tot+1,2))/(M_R(n_tot,1)-M_R(n_tot+1,1))
    error=d1_R-d1_L
    if (MOD(qn,2)==0) then
       if (error > 0.0d0) then	!Energy too high (for quantum number even)
	  go to 30
       end if
    else
       if (error < 0.0d0) then	!Energy too high (for quantum number odd)
    	  go to 30
       end if
    end if
 end do
print *,'Energy is not founded until Ei=',Ei
30	print *,'Starting shooting method.'
	print *,'Receive Ei= ',Ei,' dE= ',dE
print *,'==============================================================================='
    E1=Ei
    E2=Ei-dE
    En=(E1+E2)/2.0d0
    open(11,file='input1.txt')
    read(11,*)n_tot,E_tot,Psi0_L,Psi0_R,r0_L,r0_R,qn 
    close(11)
    do n=1,n_tot
       Psi1_R=1.0d0/10.0d0**10	!Check
       Psi1_L=1.0d0/10.0d0**10	!Check
       call numerov(n_tot,En,Psi0_L,Psi1_L,r0_L,M_L)
       r0_L=M_L(1,1)
       Psi0_L=M_L(1,2)
       call numerov(n_tot,En,Psi0_R,Psi1_R,r0_R,M_R)
       r0_R=M_R(1,1)
       Psi0_R=M_R(1,2)
       rescale=M_L(n_tot+1,2)/M_R(n_tot+1,2)
       do i=1,n_tot+2
	  M_R(i,2)=M_R(i,2)*rescale
       end do   
       d1_L=(M_L(n_tot+2,2)-M_L(n_tot+1,2))/(M_L(n_tot+2,1)-M_L(n_tot+1,1))
       d1_R=(M_R(n_tot,2)-M_R(n_tot+1,2))/(M_R(n_tot,1)-M_R(n_tot+1,1))
       error=d1_R-d1_L   
       if (abs(error) < tol) then
	  go to 40
       end if
       if (MOD(qn,2)==0) then
          if (error > 0.0d0) then	!Energy too high (for quantum number even)
	     E1=En
    	     E2=E2
	  else
             E1=E1
    	     E2=En
    	  end if
       else
          if (error < 0.0d0) then	!Energy too high (for quantum number odd)
	     E1=En
    	     E2=E2
	  else
             E1=E1
             E2=En
 	  end if
       end if
       En=(E1+E2)/2.0d0
    end do
    print *,'Energy is not founded until En=',En,'end?'
40	print*,'Starting normalization. The best En=',En
print *,'==============================================================================='
    sumL=0.0d0
    sumR=0.0d0
    do n=1,n_tot+1
	n_integL=(M_L(n,2)**2+M_L(n+1,2)**2)/2.0d0
	sumL=sumL+n_integL
        n_integR=(M_R(n,2)**2+M_R(n+1,2)**2)/2.0d0
	sumR=sumR+n_integR
    end do
    areaL=sumL*(M_L(2,1)-M_L(1,1))
    areaR=sumR*(M_R(2,1)-M_R(1,1))
    area=areaL+areaR
    norm_coef=sqrt(1.0d0/area)
    open(12,file='Rad_WF.txt')
    do n=1,n_tot+2
	M_L(n,2)=M_L(n,2)*norm_coef
        write(12,*)M_L(n,1),M_L(n,2)
    end do
    do n=1,n_tot+2
        M_R(n_tot+3-n,2)=M_R(n_tot+3-n,2)*norm_coef
        write(12,*)M_R(n_tot+3-n,1),M_R(n_tot+3-n,2)
    end do
    close(12)
print*,'Calculating the radial distribution function.'
print *,'==============================================================================='
    do n=1,n_tot+2
	M_L(n,2)=(M_L(n,2)*M_L(n,1))**2
        M_R(n,2)=(M_R(n,2)*M_R(n,1))**2
    end do
    sumL=0.0d0	!Normalize the function
    sumR=0.0d0
    do n=1,n_tot+1
      	n_integL=(M_L(n,2)+M_L(n+1,2))/2.0d0
	sumL=sumL+n_integL
        n_integR=(M_R(n,2)+M_R(n+1,2))/2.0d0
	sumR=sumR+n_integR
    end do
    areaL=sumL*(M_L(2,1)-M_L(1,1))
    areaR=sumR*(M_R(2,1)-M_R(1,1))
    area=areaL+areaR
    norm_coef=1.0d0/area
    open(13,file='Rad_distFunc.txt')
    do n=1,n_tot+2
	M_L(n,2)=M_L(n,2)*norm_coef
        write(13,*)M_L(n,1),M_L(n,2)
   end do
   do n=1,n_tot+2
        M_R(n_tot+3-n,2)=M_R(n_tot+3-n,2)*norm_coef
        write(13,*)M_R(n_tot+3-n,1),M_R(n_tot+3-n,2)
   end do
   close(13)        
print *,'Finished.'
deallocate(M_L,M_R)
end program BSwave_1D


subroutine numerov(n,E,Psi0,Psi1,r0,M_sub)
implicit none
real*8 :: ma,mb,r0,Psi0,Psi1,r_end,rmass,k2,E,re,eps
real*8 :: atr,rep,Vpot2,Gfunc2,Ffunc2,Ffunc1,r2,Upot2,Psi2
real*8 :: r1,Vpot0,Upot0,Gfunc0,Ffunc0,dr,Vpot1,Upot1,Gfunc1
real*8, dimension(n+2,2) :: M_sub
integer :: i,n,l
    open(10,file='input.txt')
    read(10,*)eps,re,ma,mb,r_end,l
    close(10)
    rmass=ma*mb/(ma+mb)
    k2=2.0d0*rmass*E
    dr=(r_end-r0)/n
    r1=r0+dr
    M_sub(1,1)=r0
    M_sub(1,2)=Psi0
    M_sub(2,1)=r1
    M_sub(2,2)=Psi1
    rep=eps*(re**12)
    atr=2*eps*(re**6)
    Vpot0=(rep/r0**12)-(atr/r0**6)
    Upot0=2.0d0*rmass*Vpot0+(l*(l+1.0d0)/r0**2)
    Gfunc0=k2-Upot0
    Ffunc0=1.0d0+Gfunc0*(dr**2/12.0d0)
    Vpot1=(rep/r1**12)-(atr/r1**6)
    Upot1=2.0d0*rmass*Vpot1+(l*(l+1.0d0)/r1**2)
    Gfunc1=k2-Upot1  
    Ffunc1=1.0d0+Gfunc1*(dr**2/12.0d0)	
    do i=1,n
     r2=r1+dr
     Vpot2=(rep/r2**12)-(atr/r2**6)
     Upot2=2.0d0*rmass*Vpot2+(l*(l+1.0d0)/r2**2)
     Gfunc2=k2-Upot2
     Ffunc2=1.0d0+Gfunc2*(dr**2/12.0d0)        
     Psi2=(Psi1*(12.0d0-10.0d0*Ffunc1)-Psi0*Ffunc0)/Ffunc2
     M_sub(i+2,1)=r2
     M_sub(i+2,2)=Psi2
	r0=r1	
	Vpot0=(rep/r0**12)-(atr/r0**6)
    	Upot0=2.0d0*rmass*Vpot0+(l*(l+1.0d0)/r0**2)
    	Gfunc0=k2-Upot0
   	Ffunc0=1.0d0+Gfunc0*(dr**2/12.0d0)
    	Psi0=Psi1
    	r1=r2
 	Vpot1=(rep/r1**12)-(atr/r1**6)
    	Upot1=2.0d0*rmass*Vpot1+(l*(l+1.0d0)/r1**2)
    	Gfunc1=k2-Upot1  
    	Ffunc1=1.0d0+Gfunc1*(dr**2/12.0d0)
    	Psi1=Psi2
    end do
return
end subroutine numerov
