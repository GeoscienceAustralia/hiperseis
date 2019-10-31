subroutine qlayer(lc, ang, n, h, valpha, vbeta, rho, qa, qb, rs, &
                w, up, wp, usv, wsv, vsh )
      implicit none
! ***********************************************
! *     crustal response for p and s waves      *
! *     in layered dissipative medium           *
! *     codedby T. Kurita & T. Mikumo           *
! *     revised by T. Shibutani                 *
! *     parallelised anf F77->F90 by I. Fomin   *
! ***********************************************
        include 'rfi_param.inc'
        integer i, k, kk, nn, m, lc, n
        real  ang
        integer,parameter :: nlmx = maxsublayer
        integer,parameter :: nb = maxdata
        integer,parameter :: nb2 = nb/2+1

      real w(nb2), h(nlmx), qa(nlmx), qb(nlmx), rs(nlmx), rho(nlmx)

      complex alpha(nlmx), ralpha(nlmx), valpha(nlmx), &
              beta(nlmx), rbeta(nlmx), vbeta(nlmx), &
              gamma(nlmx), c, rhoc2, g2,gm2, &
              en(4,4), coef1, coef2, coef3, coef4, &
              dca, dea, gca, gea, &
              up(nb2), wp(nb2), usv(nb2), wsv(nb2), &
              bm, bn, vsh(nb2), cal, cbe

      complex(8) a(4,4,nlmx), b(4,4,nlmx), aj(4,4), &
                 da, dc, de, gc, ge,  &
                 al(2,2,nlmx), bl(2,2,nlmx), ab,  &
                 p, q, cp, sp, cq, sq

      complex,parameter :: imag=(0.0,1.0)

!     *** complex velocity for dissipative medium ***

      alpha=valpha+imag*valpha/(2.0*qa)+valpha/(8.0*qa**2)
      beta=vbeta+imag*vbeta/(2.0*qb)+vbeta/(8.0*qb**2)

!     *** apparent velocity ***

      if (lc.eq.1) then
         c=alpha(n)/sin(ang)
      else
         c=beta(n)/sin(ang)
      end if

!$omp PARALLEL DO &
!$omp& default (shared) &
!$omp& private (m,cal,cbe)
      do m=1,n

        cal=(c/alpha(m))**2
        if (real(cal).ge.1.) then
           ralpha(m)=csqrt(cal-1. )
        else
           ralpha(m)=-imag*csqrt(1. -cal)
        end if

        cbe=(c/beta(m))**2
        if (real(cbe).ge.1.) then
           rbeta(m)=csqrt(cbe-1. )
        else
           rbeta(m)=-imag*csqrt(1. -cbe**2)
        end if

        gamma(m)=2.*(beta(m)/c)**2

      enddo
!$omp END PARALLEL DO

      if (lc.le.2) then

!     *** inverse matrix at the lowermost interface for P-SV ***

         en(1,1)=-2.*(beta(n)/alpha(n))**2
         en(1,2)=0.
         g2 = rho(n)*alpha(n)**2
         en(1,3)=1./g2
         en(1,4)=0.
         en(2,1)=0.
         en(2,2)=c**2*(gamma(n)-1.)/(alpha(n)**2*ralpha(n))
         en(2,3)=0.
         en(2,4)=1./(g2*ralpha(n))
         en(3,1)=(gamma(n)-1.)/(gamma(n)*rbeta(n))
         en(3,2)=0.

         g2 = rho(n)*c**2*gamma(n)
         en(3,3)=-1./(g2*rbeta(n))
         en(3,4)=0.
         en(4,1)=0.
         en(4,2)=1.
         en(4,3)=0.
         en(4,4)=1./g2

         coef1=2.*c**2/alpha(n)**2
         coef2=coef1/ralpha(n)
         coef3=2.0/gamma(n)
         coef4=coef3/rbeta(n)

!$omp PARALLEL DO &
!$omp& default (shared) &
!$omp& private (kk,nn,m,p,q,cp,sp,cq,sq,bm,a,b,aj,rhoc2) &
!$omp& private (da,dc,de,dca,dea,gc,ge,gca,gea)
      do kk=1,nb2
!       *** loop on frequency ***
        nn=n-1
        do m=1,nn
          p=w(kk)*h(m)*rs(m)/c
          q=p*rbeta(m)
          p=p*ralpha(m)
          cp=cdcos(p)
          sp=cdsin(p)
          cq=cdcos(q)
          sq=cdsin(q)
          bm=beta(m)**2*rho(m)*rbeta(m)/rs(m)
!         **  matrix at layer interfaces  **
!         ** for P-SV problem

         rhoc2=rho(m)*c**2
         a(1,1,m)=gamma(m)*cp-(gamma(m)-1.)*cq
         a(1,2,m)=imag*((gamma(m)-1.)*sp/ralpha(m)+gamma(m)*sq*rbeta(m))
         a(1,3,m)=-(cp-cq)/rhoc2
         a(1,4,m)=imag*(sp/ralpha(m)+rbeta(m)*sq)/rhoc2
         a(2,1,m)=-imag*(gamma(m)*sp*ralpha(m)+(gamma(m)-1.)*sq/rbeta(m))
         a(2,2,m)=gamma(m)*cq-(gamma(m)-1.)*cp
         a(2,3,m)=imag*(ralpha(m)*sp+sq/rbeta(m))/rhoc2
         a(3,1,m)=rhoc2*gamma(m)*(gamma(m)-1.)*(cp-cq)
         g2 = gamma(m)**2
         gm2 = (gamma(m)-1.)**2
         a(3,2,m)=imag*(rhoc2*(gm2*sp/ralpha(m)+g2*rbeta(m)*sq))
         a(4,1,m)=imag*(rhoc2*(g2*ralpha(m)*sp+gm2*sq/rbeta(m)))
         a(2,4,m)=a(1,3,m)
         a(3,3,m)=a(2,2,m)
         a(3,4,m)=a(1,2,m)
         a(4,2,m)=a(3,1,m)
         a(4,3,m)=a(2,1,m)
         a(4,4,m)=a(1,1,m)
        enddo

!       **  matrix product for P-SV case **

           m=n-1
           do while (m.ne.1)
              forall (i=1:4,k=1:4) b(i,k,m-1)=sum(a(i,1:4,m)*a(1:4,k,m-1))
              a(1:4,1:4,m-1) = b(1:4,1:4,m-1)
              m=m-1
           enddo
           forall (i=1:4,k=1:4) aj(i,k)=sum(en(i,1:4)*a(1:4,k,1))
           da= (aj(1,1)-aj(2,1))*(aj(3,2)-aj(4,2)) &
              -(aj(1,2)-aj(2,2))*(aj(3,1)-aj(4,1))
        if (lc.eq.1) then

!       **  up horizontal/total, wp vertical/total of P  **

           dc=aj(4,2)-aj(3,2)
           de=aj(4,1)-aj(3,1)
           dca=CMPLX(dc/da)
           dea=CMPLX(de/da)
           up(kk)=dca*coef1*alpha(n)/c
           wp(kk)=dea*coef2*alpha(n)*ralpha(n)/c

        else if (lc.eq.2) then

!     **  usv horizontal/total, wsv vertical/total of SV  **

           gc=aj(1,2)-aj(2,2)
           ge=aj(2,1)-aj(1,1)
           gca=CMPLX(gc/da)
           gea=CMPLX(ge/da)
           usv(kk)=gca*coef4*beta(n)*rbeta(n)/c
           wsv(kk)=gea*coef3*beta(n)/c

        end if

      enddo
!$omp END PARALLEL DO
      else !lc > 2
!$omp PARALLEL DO &
!$omp& default (shared) &
!$omp& private (kk,nn,m,p,q,cq,sq,bm,al,bl,bn,ab)
      do kk=1,nb2
!        *** loop on frequency ***
        nn=n-1
        do m=1,nn
          p=w(kk)*h(m)*rs(m)/c
          q=p*rbeta(m)
          p=p*ralpha(m)
          cq=cdcos(q)
          sq=cdsin(q)
          bm=beta(m)**2*rho(m)*rbeta(m)/rs(m)

!         **  matrix at layer interfaces for SH problem **
          al(1,1,m)=cq
          al(2,2,m)=al(1,1,m)
          al(1,2,m)=imag*sq/bm
          al(2,1,m)=imag*sq*bm
        enddo

!       **  matrix product for SH waves  **

         m=n-1
         do while (m.ne.1)
            forall(i=1:2,k=1:2) bl(i,k,m-1)=sum(al(i,1:2,m)*al(1:2,k,m-1))
            al(1:2,1:2,m-1)=bl(1:2,1:2,m-1)
            m=m-1
         enddo
         bn=beta(n)**2*rho(n)*rbeta(n)/rs(n)
         ab=al(2,1,1)+bn*al(1,1,1)

!       **  vsh horizontal/total of SH  **

         vsh(kk)=2.*cmplx(bn/ab)

      enddo
!$omp END PARALLEL DO
      endif

      end subroutine qlayer
