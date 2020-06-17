      subroutine theo(
     &     n, fbeta, h, vps, qa, qb, fs, din, a0, c0, t0,
     &     nd, rx )
      implicit none
******************************************
*     crustal response for layered dissipative medium       *
*     coded by t. shibutani    ( jan. 1990 )                *
******************************************
c
c     **  din(deg),rin(rad) = incident angle of p wave
c     **  vpvs = vp/vs
c     **  qa,qb = p and s attenuations
c     **  n = the number of layers
c     **  fbeta(m) = s velocity in the m-th layer
c     **  h(m) = thickness of the m-th layer
c     **  c0 = parameter of a minimum amplitude level
c     **  a0 = parameter of a gaussian high-cut filter
c     **  rx(1-nb) = synthetic response
c
c
      include 'rfi_param.inc'

      integer n, m, nlmx, nb, nb2, nd, i, k
      real fs, din, rad, a0, c0, t0, fk, fr, pi, pi2, frini, frint
      real rin, rvmax, va, wmax
      parameter     ( nlmx = maxsublayer)
      parameter     ( nb = freq_resolution )
      parameter     ( nb2 = nb/2+1 )

c     parameter     ( nlmx = 100)
c     parameter     ( nb = 1024 )
c     parameter     ( nb2 = nb/2+1 )

c       !!! remarks: if you change the above parameters, you should
c       !!!          also change the corresponding parameters in
c       !!!          the main program and the subroutine "qlayer".
c
c
      real w(nb2), wa(nb2), fai(nb2), gau(nb2),
     &     data(2*nb),
     &     fbeta(nlmx), h(nlmx), vps(nlmx), rho(nlmx),
     &     qa(nlmx), qb(nlmx), rs(nlmx),
     &     rx(maxdata)
c
      complex up(nb2), wp(nb2), usv(nb2), wsv(nb2), vsh(nb2),
     &        cc(nb),
     &        valpha(nlmx), vbeta(nlmx),
     &        yi, cs
c
c     *** parameters ***
      yi=(0.,1.)
      pi=3.141592653
      pi2=pi*2
      rad=180./pi
      frini=0.
      frint=fs/nb
c

c      write(*,*)'teo1',nd

      do 110 k=1,nb2
        fk=k-1
        fr=frini+frint*fk
        w(k)=pi2*fr
  110 continue
c
c     *** crustal response for layered q structure ***
c     **  p velocity,s velocity,density,thickness  **
c     **  p and s attenuations, lateral anisotropy  **
c
c	do i=1,n
c	  write(6,'(i3,5f10.3)') i,h(i),fbeta(i),vps(i),qa(i),qb(i)
c	end do
c      write(6,*) 'nd=',nd
c
      do 120 m=1,n
        vbeta(m)=cmplx(fbeta(m))
        va=vps(m)*fbeta(m)
        valpha(m)=cmplx(va)
        rho(m)=2.35+0.036*(va-3.0)**2
        rs(m)=1.0
  120 continue

c      write(*,*)'teo2'

c
c     **  computation for the case of incident p waves  **
      rin=pi*din/180.
      call qlayer(
     &     1, rin, n, h, valpha, vbeta, rho, qa, qb, rs,
     &     w, up, wp, usv, wsv, vsh )
c      write(*,*)'teo3'

c
c     * a minimum allowable amplitude level *
      wmax=0.
      do 40 i=1,nb2
        wa(i)=real(wp(i)*conjg(wp(i)))
        wmax=amax1(wmax,wa(i))
   40 continue
      do 50 i=1,nb2
        fai(i)=amax1(wa(i),c0*wmax)
   50 continue
c
c     * a gaussian function to exclude high-frequency noise *
      do 60 i=1,nb2
        gau(i)=exp(-(w(i)/(2*a0))**2)
   60 continue
c
c     * vertical receiver function for absolute amplitude *
c
      do 260 i=1,nb2
        cs=yi*w(i)*t0
        cc(i)=wp(i)*conjg(wp(i))*gau(i)*cexp(cs)/fai(i)
        if (i.ge.2) cc(nb-i+2)=conjg(cc(i))
  260 continue
c	write(*,*)'teo4'
c     * transform back into the time domain *
c
      do 270 i=1,nb
        data(2*i-1)=real(cc(i))
        data(2*i)=aimag(cc(i))
  270 continue
      call four1(data,nb,-1)
c	write(*,*)'teo5'

c
c     * max amplitude of vertical reciever function *
c
      rvmax=0.
      do 510 i=1,nd
        rvmax=amax1(rvmax,abs(data(2*i-1)/nb))

  510 continue
c
c     * the receiver function in the frequency domain *
c
      do 160 i=1,nb2
        cs=yi*w(i)*t0
        cc(i)=conjg(up(i))*wp(i)*gau(i)*cexp(cs)/fai(i)
        if (i.ge.2) cc(nb-i+2)=conjg(cc(i))
  160 continue

c
c     * transform back into the time domain *
c
      do 170 i=1,nb
        data(2*i-1)=real(cc(i))
        data(2*i)=aimag(cc(i))
  170 continue
      call four1(data,nb,-1)
c	write(*,*)'teo6'
c
c     **  radial receiver function ***
c
      do 520 i=1,nd
        rx(i)=data(2*i-1)/nb/rvmax
c        write(*,*)'teo6b',i,nd
  520 continue
c      write(6,*) 'rx(71:80)'
c      write(6,*) (rx(i),i=71,80)
c
c	write(*,*)'teo7'
      return
      end
