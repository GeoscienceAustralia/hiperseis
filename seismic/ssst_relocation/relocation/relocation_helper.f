!file: relocation_helper.f

      subroutine relocation(elon0,ecolat0,edep0,slon0,scolat0,selev0,
     &                      iph0,phases_temp,wt0,tt,term,dlon,dcolat,
     &                      ddep,niter,frac,inorm,ttarr,dtddarr,xarr1,
     &                      darr1,xarr2,darr2,tau0,tau1,tau2,npick_use,
     &                      elonbest,ecolatbest,edepbest,tbest,resid,
     &                      qual,
     &                      npick,nx1,nd1,nx2,nd2,nph)

      implicit none

! Dependent on other variables

      integer npick
Cf2py intent(in) npick
      integer nx1
Cf2py intent(in) nx1
      integer nd1
Cf2py intent(in) nd1
      integer nx2
Cf2py intent(in) nx2
      integer nd2
Cf2py intent(in) nd2
      integer nph
Cf2py intent(in) nph

! Input variables

      real elon0
Cf2py intent(in) elon0
      real ecolat0
Cf2py intent(in) ecolat0
      real edep0
Cf2py intent(in) edep0
      real slon0(npick)
Cf2py intent(in) slon0
      real scolat0(npick)
Cf2py intent(in) scolat0
      real selev0(npick)
Cf2py intent(in) selev0
      integer iph0(npick)
Cf2py intent(in) iph0
      character phases_temp(npick,8)
Cf2py intent(in) phases_temp
      integer wt0(npick)
Cf2py intent(in) wt0
      real tt(npick)
Cf2py intent(in) tt
      real term(npick)
Cf2py intent(in) term
      real dlon
Cf2py intent(in) dlon
      real dcolat
Cf2py intent(in) dcolat
      real ddep
Cf2py intent(in) ddep
      integer niter
Cf2py intent(in) niter
      real frac
Cf2py intent(in) frac
      integer inorm
Cf2py intent(in) inorm
      real ttarr(nx1,nd1,nph)
Cf2py intent(in) ttarr
      real dtddarr(nx1,nd1,nph)
Cf2py intent(in) dtddarr
      real tau0(nx2,nd2,nph)
Cf2py intent(in) tau0
      real tau1(nx2,nd2,nph)
Cf2py intent(in) tau1
      real tau2(nx2,nd2,nph)
Cf2py intent(in) tau2
      real xarr1(nx1,nph)
Cf2py intent(in) xarr1
      real darr1(nd1,nph)
Cf2py intent(in) darr1
      real xarr2(nx2,nph)
Cf2py intent(in) xarr2
      real darr2(nd2,nph)
Cf2py intent(in) darr2
      integer npick_use
Cf2py intent(in) npick_use

! Output variables

      real elonbest
Cf2py intent(out) elonbest
      real ecolatbest
Cf2py intent(out) ecolatbest
      real edepbest
Cf2py intent(out) edepbest
      real tbest
Cf2py intent(out) tbest
      real resid(npick)
Cf2py intent(out) resid
      integer qual
Cf2py intent(out) qual

! Parameters

      real degkm
      parameter (degkm=111.19493)
      real degrad
      parameter (degrad=57.29578)

! Other variables

      real pi
      real elon
      real ecolat
      real edep
      real fitbest
      integer i
      integer j
      character (len=8) phases(npick)
      integer it
      integer ix
      integer iy
      integer iz
      integer npick_used
      integer ip
      real slon
      real scolat
      real selev
      integer iph
      character*8 phase
      integer wt
      real ecdist
      real angdist
      real azim
      integer iflag
      real ttime
      real ellipcorr
      real elevcorr
      real resid2(npick)
      real fit
      real residval

! Set up variables

      pi = 4*atan(1.0)
      ecolatbest=ecolat0
      elonbest=elon0
      edepbest=edep0
      ecolat=ecolatbest
      elon=elonbest
      edep=edepbest
      tbest=0.0
      fitbest=9.e20

      do i=1,npick
         do j=1,8
            phases(i)(j:j) = phases_temp(i,j)
         end do
      end do

! relocation algorithm

      do it=1,niter
         elon0=elonbest
         ecolat0=ecolatbest
         edep0=edepbest
!         print *, 'Iteration ',it
!         print *, elon0,ecolat0,edep0,tbest,fitbest
         do 1 ix=-1,1
            elon=elon0+dlon*float(ix)
            if(elon.lt.-180.0) then
               elon=elon+360.0
            else if(elon.gt.180.0) then
               elon=elon-360.0
            end if
            elon=elon/degrad
            do 2 iy=-1,1
               ecolat=ecolat0+dcolat*float(iy)
               if((ecolat.lt.0).or.(ecolat.gt.180.0)) then
                  go to 2
               end if
               ecolat=ecolat/degrad
               do 3 iz=-1,1
                  edep=edep0+ddep*float(iz)
                  if(edep.lt.0) then
                     go to 3
                  end if
                  npick_used=0
                  ip=1
!                  print *, ix,iy,iz
                  do 4 while ((ip.le.npick).and.
     &                        (npick_used.lt.npick_use))
                     slon=slon0(ip)/degrad
                     scolat=scolat0(ip)/degrad
                     selev=selev0(ip)
                     iph=iph0(ip)
                     phase=phases(ip)
                     wt=wt0(ip)
                     ecdist=acos(cos(ecolat)*cos(scolat)+
     &                           sin(ecolat)*sin(scolat)*
     &                           cos(elon-slon))
                     angdist=ecdist*degrad

                     call azimuth(elon,pi/2.0-ecolat,slon,
     &                            pi/2.0-scolat,azim)
                     call GET_TT(ttarr(:,:,iph),dtddarr(:,:,iph),
     &                           xarr1(:,iph),darr1(:,iph),angdist,
     &                           edep,selev,wt,ttime,elevcorr,iflag,
     &                           nx1,nd1)
                     call ellip_corr(tau0(:,:,iph),tau1(:,:,iph),
     &                               tau2(:,:,iph),xarr2(:,iph),
     &                               darr2(:,iph),angdist,edep,azim,
     &                               ecolat,ellipcorr,nx2,nd2)

                     if(iflag.eq.-1.or.ttime.eq.999.0) then
                        resid(ip)=-999.0
                        ip=ip+1
                        go to 4
                     end if

                     if(abs(elevcorr).ge.1.0) then
                        elevcorr=0.0
                     end if
                     npick_used=npick_used+1
                     ttime=ttime+elevcorr+ellipcorr
                     resid(ip)=tt(ip)-term(ip)-ttime
                     resid2(npick_used)=resid(ip)
                     ip=ip+1
4                 continue

                  if(npick_used.le.0) then
                     fit=10.e20
                  else
                     fit=0.0
                     if(inorm.eq.1) then
                        call median(resid2,npick_used,residval)
                        do i=1,npick_used
                           fit=fit+abs(resid2(i)-residval)
                        end do
                        fit=fit/float(npick_used)
                     else if(inorm.eq.2) then
                        call mean(resid2,npick_used,residval)
                        do i=1,npick_used
                           fit=fit+(resid2(i)-residval)**2
                        end do
                        fit=sqrt(fit/npick_used)
                     end if
                  end if

!                  print *, elon*degrad,ecolat*degrad,edep,residval,
!     &                     fit
                  if(fit.lt.fitbest) then
                     fitbest=fit
                     elonbest=elon*degrad
                     ecolatbest=ecolat*degrad
                     edepbest=edep
                     tbest=residval
!                     print *, 'New hypocentre'
                  end if
3              continue
2           continue
1        continue

         dlon=dlon*frac
         dcolat=dcolat*frac
         ddep=ddep*frac
      end do

! residual calculation

      npick_used=0
      elon=elonbest/degrad
      ecolat=ecolatbest/degrad
      edep=edepbest
      do 5 ip=1,npick
         slon=slon0(ip)/degrad
         scolat=scolat0(ip)/degrad
         selev=selev0(ip)
         iph=iph0(ip)
         phase=phases(ip)
         wt=wt0(ip)
         ecdist=acos(cos(ecolat)*cos(scolat)+
     &               sin(ecolat)*sin(scolat)*cos(elon-slon))
         angdist=ecdist*degrad

         call azimuth(elon,pi/2.0-ecolat,slon,pi/2.0-scolat,azim)
         call get_tt(ttarr(:,:,iph),dtddarr(:,:,iph),xarr1(:,iph),
     &               darr1(:,iph),angdist,edep,selev,wt,ttime,elevcorr,
     &               iflag,nx1,nd1)
         call ellip_corr(tau0(:,:,iph),tau1(:,:,iph),tau2(:,:,iph),
     &                   xarr2(:,iph),darr2(:,iph),angdist,edep,azim,
     &                   ecolat,ellipcorr,nx2,nd2)

         if(iflag.eq.-1.or.ttime.eq.999.0) then
            resid(ip)=-999.0
            go to 5
         end if

         if(abs(elevcorr).ge.1.0) then
            elevcorr=0.0
         end if
         npick_used=npick_used+1
         ttime=ttime+elevcorr+ellipcorr
         resid(ip)=tt(ip)-term(ip)-ttime-tbest
5     continue

      if(npick_used.le.3) then
         qual=-1
      end if

      return
      end

!----------------------------------------------------------------------
     
      subroutine get_tt(t,dtdd,x,d,del,dep,selev,wt,tt,elevcorr,iflag,
     &                  nx0,nd0)

      implicit none

! Dependent on other variables

      integer nd0
Cf2py intent(in) nd0
      integer nx0
Cf2py intent(in) nx0

! Input variables

      real t(nx0,nd0)
Cf2py intent(in) t
      real dtdd(nx0,nd0)
Cf2py intent(in) dtdd
      real x(nx0)
Cf2py intent(in) x
      real d(nd0)
Cf2py intent(in) d
      real del
Cf2py intent(in) del
      real dep
Cf2py intent(in) dep
      real selev
Cf2py intent(in) selev
      integer wt
Cf2py intent(in) wt

! Output variables

      real tt
Cf2py intent(out) tt
      real elevcorr
Cf2py intent(out) elevcorr
      integer iflag
Cf2py intent(out) iflag

! Other variables

      integer id0
      integer id1
      integer id2
      integer id3
      integer ix0
      integer ix1
      integer ix2
      integer ix3
      integer nd
      integer nx
      integer flag(10)
      integer iflag1
      integer iflag2

      real dfrac
      real ttemp(4)
      real t1
      real t2
      real surfvel
      real dtddtemp(4)
      real dtdd1
      real dtdd2
      real dtdd_val
      real xfrac
      real psurfvel
      real ssurfvel
      real degkm

! Parameters

      parameter (psurfvel=5.80)
      parameter (ssurfvel=3.46)
      parameter (degkm=111.19493)

! Check if P/S wave

      if(wt.eq.1) then
         surfvel=psurfvel
      else
         surfvel=ssurfvel
      end if

! Check if outside coordinate range
! Note that id1 < 1, id1 > nd, ix1 < 1, and ix1 > nx are all caught by
! the nearest_indices subroutine and the flag will return -1 if any of
! these cases are true.

      nx = maxloc(x,1)
      nd = maxloc(d,1)

      call nearest_indices(x(1:nx),nx,del,ix1,ix2,iflag1)
      call nearest_indices(d(1:nd),nd,dep,id1,id2,iflag2)

      id0=id1-1
      id3=id2+1
      ix0=ix1-1
      ix3=ix2+1

      if((iflag1.eq.-1).or.(iflag2.eq.-1)) then
         go to 3
      else if(((id1.le.1).or.(id2.ge.nd)).or.
     &        ((ix1.le.1).or.(ix2.ge.nx))) then
         go to 2
      else
         go to 1
      end if

! Interpolate to find travel time
! Spline interpolation

1     if((isnan(t(ix0,id0)).or.isnan(t(ix0,id3))).or.
     &   (isnan(t(ix3,id0)).or.isnan(t(ix3,id3)))) then
         go to 2
      end if

      iflag=0

      call spline(x(ix0:ix3),t(ix0:ix3,id0),3,del,ttemp(1),
     &            flag(1))
      call spline(x(ix0:ix3),t(ix0:ix3,id1),3,del,ttemp(2),
     &            flag(2))
      call spline(x(ix0:ix3),t(ix0:ix3,id2),3,del,ttemp(3),
     &            flag(3))
      call spline(x(ix0:ix3),t(ix0:ix3,id3),3,del,ttemp(4),
     &            flag(4))
      call spline(d(id0:id3),ttemp,3,dep,tt,flag(5))

      call spline(x(ix0:ix3),dtdd(ix0:ix3,id0),3,del,dtddtemp(1),
     &            flag(6))
      call spline(x(ix0:ix3),dtdd(ix0:ix3,id1),3,del,dtddtemp(2),
     &            flag(7))
      call spline(x(ix0:ix3),dtdd(ix0:ix3,id2),3,del,dtddtemp(3),
     &            flag(8))
      call spline(x(ix0:ix3),dtdd(ix0:ix3,id3),3,del,dtddtemp(4),
     &            flag(9))
      call spline(d(id0:id3),dtddtemp,3,dep,dtdd_val,flag(10))

      if(sum(flag).lt.0) then
         go to 2
      end if

      elevcorr=(surfvel*dtdd_val/degkm)**2
      if (elevcorr.gt.1.0) then
         elevcorr=1.0/elevcorr
      end if
      elevcorr=sqrt(1.0-elevcorr)
      elevcorr=elevcorr*selev/(1000*surfvel)

      return

! Bilinear interpolation

2     if((isnan(t(ix1,id1)).or.isnan(t(ix1,id2))).or.
     &   (isnan(t(ix2,id1)).or.isnan(t(ix2,id2)))) then
         go to 3
      end if

      iflag=0
      xfrac=(del-x(ix1))/(x(ix2)-x(ix1))
      dfrac=(dep-d(id1))/(d(id2)-d(id1))

      t1=t(ix1,id1)+xfrac*(t(ix2,id1)-t(ix1,id1))
      t2=t(ix1,id2)+xfrac*(t(ix2,id2)-t(ix1,id2))
      tt=t1+dfrac*(t2-t1)

      dtdd1=dtdd(ix1,id1)+xfrac*(dtdd(ix2,id1)-dtdd(ix1,id1))
      dtdd2=dtdd(ix1,id2)+xfrac*(dtdd(ix2,id2)-dtdd(ix1,id2))
      dtdd_val=dtdd1+dfrac*(dtdd2-dtdd1)
      elevcorr=(surfvel*dtdd_val/degkm)**2
      if (elevcorr.gt.1.0) then
         elevcorr=1.0/elevcorr
      end if
      elevcorr=sqrt(1.0-elevcorr)
      elevcorr=elevcorr*selev/(1000*surfvel)

      return

! No travel time

3     iflag=-1
      tt=999.0
      elevcorr=999.0

      return
      end

!----------------------------------------------------------------------
     
      subroutine ellip_corr(tau0,tau1,tau2,x,d,del,dep,azim,colat,
     &                      ellipcorr,nx0,nd0)

      implicit none

! Dependent on other variables

      integer nd0
Cf2py intent(in) nd0
      integer nx0
Cf2py intent(in) nx0

! Input variables

      real tau0(nx0,nd0)
Cf2py intent(in) tau0
      real tau1(nx0,nd0)
Cf2py intent(in) tau1
      real tau2(nx0,nd0)
Cf2py intent(in) tau2
      real x(nx0)
Cf2py intent(in) x
      real d(nd0)
Cf2py intent(in) d
      real del
Cf2py intent(in) del
      real dep
Cf2py intent(in) dep
      real azim
Cf2py intent(in) azim
      real colat
Cf2py intent(in) colat

! Output variables

      real ellipcorr
Cf2py intent(out) ellipcorr

! Other variables

      integer id0
      integer id1
      integer id2
      integer id3
      integer ix0
      integer ix1
      integer ix2
      integer ix3
      integer nd
      integer nx
      integer flag(15)
      integer iflag1
      integer iflag2

      real dfrac
      real xfrac
      real temp(4)
      real v1
      real v2
      real sc0
      real sc1
      real sc2
      real tau0val
      real tau1val
      real tau2val

! Check if outside coordinate range
! Note that id1 < 1, id1 > nd, ix1 < 1, and ix1 > nx are all caught by
! the nearest_indices subroutine and the flag will return -1 if any of
! these cases are true.

      nx = maxloc(x,1)
      nd = maxloc(d,1)

      call nearest_indices(x(1:nx),nx,del,ix1,ix2,iflag1)
      call nearest_indices(d(1:nd),nd,dep,id1,id2,iflag2)

      id0=id1-1
      id3=id2+1
      ix0=ix1-1
      ix3=ix2+1

      sc0 = (1.0+3.0*cos(2.0*colat))/4.0
      sc1 = sqrt(3.0)*sin(2.0*colat)/2.0
      sc2 = sqrt(3.0)*sin(colat)*sin(colat)/2.0

      if((iflag1.eq.-1).or.(iflag2.eq.-1)) then
         go to 3
      else if(((id1.le.1).or.(id2.ge.nd)).or.
     &        ((ix1.le.1).or.(ix2.ge.nx))) then
         go to 2
      else
         go to 1
      end if

! Interpolate to find time correction
! Spline interpolation

1     call spline(x(ix0:ix3),tau0(ix0:ix3,id0),3,del,temp(1),
     &            flag(1))
      call spline(x(ix0:ix3),tau0(ix0:ix3,id1),3,del,temp(2),
     &            flag(2))
      call spline(x(ix0:ix3),tau0(ix0:ix3,id2),3,del,temp(3),
     &            flag(3))
      call spline(x(ix0:ix3),tau0(ix0:ix3,id3),3,del,temp(4),
     &            flag(4))
      call spline(d(id0:id3),temp,3,dep,tau0val,flag(5))

      call spline(x(ix0:ix3),tau1(ix0:ix3,id0),3,del,temp(1),
     &            flag(6))
      call spline(x(ix0:ix3),tau1(ix0:ix3,id1),3,del,temp(2),
     &            flag(7))
      call spline(x(ix0:ix3),tau1(ix0:ix3,id2),3,del,temp(3),
     &            flag(8))
      call spline(x(ix0:ix3),tau1(ix0:ix3,id3),3,del,temp(4),
     &            flag(9))
      call spline(d(id0:id3),temp,3,dep,tau1val,flag(10))

      call spline(x(ix0:ix3),tau2(ix0:ix3,id0),3,del,temp(1),
     &            flag(11))
      call spline(x(ix0:ix3),tau2(ix0:ix3,id1),3,del,temp(2),
     &            flag(12))
      call spline(x(ix0:ix3),tau2(ix0:ix3,id2),3,del,temp(3),
     &            flag(13))
      call spline(x(ix0:ix3),tau2(ix0:ix3,id3),3,del,temp(4),
     &            flag(14))
      call spline(d(id0:id3),temp,3,dep,tau2val,flag(15))

      if(sum(flag).lt.0) then
         go to 2
      end if

      ellipcorr=sc0*tau0val+sc1*cos(azim)*tau1val+
     &          sc2*cos(2.0*azim)*tau2val

      return

! Bilinear interpolation

2     xfrac=(del-x(ix1))/(x(ix2)-x(ix1))
      dfrac=(dep-d(id1))/(d(id2)-d(id1))

      v1=tau0(ix1,id1)+xfrac*(tau0(ix2,id1)-tau0(ix1,id1))
      v2=tau0(ix1,id2)+xfrac*(tau0(ix2,id2)-tau0(ix1,id2))
      tau0val=v1+dfrac*(v2-v1)

      v1=tau1(ix1,id1)+xfrac*(tau1(ix2,id1)-tau1(ix1,id1))
      v2=tau1(ix1,id2)+xfrac*(tau1(ix2,id2)-tau1(ix1,id2))
      tau1val=v1+dfrac*(v2-v1)

      v1=tau2(ix1,id1)+xfrac*(tau2(ix2,id1)-tau2(ix1,id1))
      v2=tau2(ix1,id2)+xfrac*(tau2(ix2,id2)-tau2(ix1,id2))
      tau2val=v1+dfrac*(v2-v1)

      ellipcorr=sc0*tau0val+sc1*cos(azim)*tau1val+
     &          sc2*cos(2.0*azim)*tau2val

      return

! No time correction

3     ellipcorr=0.0

      return
      end

!----------------------------------------------------------------------

      subroutine azimuth(elon,elat,slon,slat,azim)

      implicit none

      real elon
      real elat
      real slon
      real slat
      real azim
Cf2py intent(in) elon
Cf2py intent(in) elat
Cf2py intent(in) slon
Cf2py intent(in) slat
Cf2py intent(out) azim

      real pi
      pi=4*atan(1.0)

      azim=atan(sin(slon-elon)/
     &          (cos(elat)*tan(slat)-sin(elat)*cos(slon-elon)))

      if (slat.lt.elat) then
         azim=azim+pi
      end if
      if ((slat.eq.elat).and.(elat.lt.0.0)) then
         azim=azim+pi
      end if
      if (azim.lt.0.0) then
         azim=azim+2*pi
      else if (azim.gt.2.0*pi) then
         azim=azim-2*pi
      end if

      return
      end

!----------------------------------------------------------------------

      subroutine mean(x,n,xmean)

      implicit none

      integer n
Cf2py intent(in) n

      real x(n)
Cf2py intent(in) x

      real xmean
Cf2py intent(out) xmean

      integer i
      real s

      if (n.eq.0) then
         xmean=0.
         return
      end if

      s=0.0
      do i=1,n
         s=s+x(i)
      end do
      xmean=s/float(n)
      return
      end

!----------------------------------------------------------------------

      subroutine median(x,n,xmed)

      implicit none

      integer n
Cf2py intent(in) n

      real x(n)
Cf2py intent(in) x

      real xmed
Cf2py intent(out) xmed

      if(n.eq.0) then
         xmed=0.0
         return
      else if(n.eq.1) then
         xmed=x(1)
      else
         call quicksort(x,n)
         if(mod(n,2).eq.0) then
            xmed=(x(n/2)+x(n/2+1))/2.0
         else
            xmed=x((n+1)/2)
         end if
      end if
      return
      end
      
!----------------------------------------------------------------------

      recursive subroutine quicksort(x,n)

      implicit none

      integer n
Cf2py intent(in) n
      real x(n)
Cf2py intent(in, out) x

      integer ne
      real xe(n)
      integer nl
      real xl(n)
      integer nh
      real xh(n)
      
      integer i

      if(n.gt.1) then
         ne=0
         nl=0
         nh=0
         do i=1,n
            if(x(i).lt.x(1)) then
               nl=nl+1
               xl(nl)=x(i)
            else if(x(i).gt.x(1)) then
               nh=nh+1
               xh(nh)=x(i)
            else
               ne=ne+1
               xe(ne)=x(i)
            end if
         end do

         call quicksort(xl(:nl),nl)
         call quicksort(xh(:nl),nh)

         do i=1,nl
            x(i)=xl(i)
         end do
         do i=1,ne
            x(i+nl)=xe(i)
         end do
         do i=1,nh
            x(i+nl+ne)=xh(i)
         end do
      end if
      return
      end

!----------------------------------------------------------------------

      subroutine spline(t,a,n,x,y,flag)

      implicit none

      integer n
Cf2py intent(in) n
      real t(n+1)
Cf2py intent(in) t
      real a(n+1)
Cf2py intent(in) a
      real x
Cf2py intent(in) x
      real y
Cf2py intent(out) y
      integer flag
Cf2py intent(out) flag

      real b(n)
      real c(n+1)
      real d(n)
      real h(n)
      real k(n)
      real l(n+1)
      real m(n+1)
      real z(n+1)
      integer i
      integer ilo
      integer ihi

      do i=1,n
         h(i)=t(i+1)-t(i)
      end do

      k(1)=0.0
      do i=2,n
         k(i)=(3.0*(a(i+1)-a(i))/h(i))-(3.0*(a(i)-a(i-1))/h(i-1))
      end do

      l(1)=1.0
      m(1)=0.0
      z(1)=0.0

      do i=2,n
         l(i)=2.0*(t(i+1)-t(i-1))-h(i-1)*m(i-1)
         m(i)=h(i)/l(i)
         z(i)=(k(i)-h(i-1)*z(i-1))/l(i)
      end do

      l(n+1)=1.0
      c(n+1)=0.0
      z(n+1)=0.0

      do i=n,1,-1
         c(i)=z(i)-m(i)*c(i+1)
         b(i)=((a(i+1)-a(i))/h(i))-(h(i)*(c(i+1)+2.0*c(i))/3.0)
         d(i)=(c(i+1)-c(i))/(3.0*h(i))
      end do

      call nearest_indices(t,n+1,x,ilo,ihi,flag)

      if(flag.eq.-1) then
         y=0.0
         return
      end if

      y=a(ilo)+b(ilo)*(x-t(ilo))+c(ilo)*(x-t(ilo))**2+
     &  d(ilo)*(x-t(ilo))**3

      return
      end

!----------------------------------------------------------------------

      subroutine nearest_indices(x,n,y,ilo,ihi,flag)

      implicit none

      integer n
Cf2py intent(in) n
      real x(n)
Cf2py intent(in) x
      real y
Cf2py intent(in) y
      integer ilo
Cf2py intent(out) ilo
      integer ihi
Cf2py intent(out) ihi
      integer flag
Cf2py intent(out) flag

      integer i

      flag=1
      ilo=1
      ihi=n

      if((y.lt.x(ilo)).or.(y.gt.x(ihi))) then
         flag=-1
      end if

      do while(flag.eq.1)
         if((ihi-ilo).le.1) then
            flag=0
         else
            if(mod(ihi+1-ilo,2).eq.0) then
               i=(ihi-1+ilo)/2
            else
               i=(ihi+ilo)/2
            end if
            if(y.gt.x(i)) then
               ilo=i
            else
               ihi=i
            end if
         end if
      end do
      return
      end

!end file: relocation_helper.f