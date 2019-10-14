program test
    implicit none
    include 'mpif.h'
    real , EXTERNAL    ::    gasdev,ran3
    real log, sqrt
    integer v
    do
        read(*, *) v
        if (v == 0.0) then
            exit
        endif
        write(*, *) v, gasdev(v), ran3(v)
    enddo
end program test



!-------------------------------------------------------------------
!                       
!   Numerical Recipes random number generator for 
!       a Gaussian distribution
!
! ----------------------------------------------------------------------------



FUNCTION GASDEV(idum)
implicit none
!     ..Arguments..
integer          idum
real GASDEV

!     ..Local..
real v1,v2,r,fac
real ran3

if (idum.lt.0) iset=0
10   v1=2*ran3(idum)-1
v2=2*ran3(idum)-1
r=v1**2+v2**2
if(r.ge.1.or.r.eq.0) GOTO 10
fac=sqrt(-2*log(r)/r)
GASDEV=v2*fac

RETURN
END


!-------------------------------------------------------------------
!                       
!   Numerical Recipes random number generator
!
! ----------------------------------------------------------------------------

FUNCTION ran3(idum)
implicit none
INTEGER idum
INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
! write(*,*)' idum ',idum
if(idum.lt.0.or.iff.eq.0)then
    iff=1
    mj=MSEED-iabs(idum)
    mj=mod(mj,MBIG)
    ma(55)=mj
    mk=1
    do 11 i=1,54
    ii=mod(21*i,55)
    ma(ii)=mk
    mk=mj-mk
    if(mk.lt.MZ)mk=mk+MBIG
    mj=ma(ii)
!  write(*,*)' idum av',idum
11      continue
    do 13 k=1,4
    do 12 i=1,55
    ma(i)=ma(i)-ma(1+mod(i+30,55))
    if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
! write(*,*)' idum ap',idum
    inext=0
    inextp=31
    idum=1
endif
! write(*,*)' idum app ',idum
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
!  write(*,*)' idum ',idum
    
return
END

!the code is finished. 
