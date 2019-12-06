      subroutine qlayer( 
     &     lc, ang, n, h, valpha, vbeta, rho, qa, qb, rs, 
     &     w, up, wp, usv, wsv, vsh )
      use, intrinsic :: iso_fortran_env
      implicit none
***********************************************                         
*     crustal response for p and s waves      *                         
*     in layered dissipative medium           *                         
*     codedby T. Kurita & T. Mikumo           *                         
*     revised by T. Shibutani                 *                         
***********************************************                         
c                                                                       
c
        include 'rfi_param.inc'
        integer nlmx, nb, nb2, i, j, k, kk, mm, nn, m, lc, n
        real  ang, cal, cbe
      parameter     ( nlmx = maxsublayer )
      parameter     ( nb  = maxdata )
      parameter     ( nb2 = nb/2+1 )
c
c
      real w(nb2), h(nlmx), qa(nlmx), qb(nlmx), rs(nlmx), rho(nlmx)
c
      complex alpha(nlmx), ralpha(nlmx), valpha(nlmx),
     &        beta(nlmx), rbeta(nlmx), vbeta(nlmx),
     &        gamma(nlmx), vn, c, imag, 
     &        en(4,4), coef1, coef2, coef3, coef4,
     &        dca, dea, gca, gea,
     &        up(nb2), wp(nb2), usv(nb2), wsv(nb2),
     &        bm, bn, rvs, vsh(nb2)
c                                                                       
      complex(kind=real64) a(4,4,nlmx), b(4,4,nlmx), aj(4,4),
     &           s, da, dc, de, gc, ge, 
     &           al(2,2,nlmx), bl(2,2,nlmx), sl, ab, 
     &           p(nlmx), q(nlmx), cp, sp, cq, sq
c
      parameter ( imag=(0.0,1.0) )                                      
c                                                                       
c     *** complex velocity for dissipative medium ***
c
      do 1 m=1,n                                                        
        alpha(m)=valpha(m)+imag*valpha(m)/(2.0*qa(m))                   
     &          +valpha(m)/(8.0*qa(m)**2)                               
        beta(m)=vbeta(m)+imag*vbeta(m)/(2.0*qb(m))                      
     &         +vbeta(m)/(8.0*qb(m)**2)                                 
    1 continue                                                          
c
      if (lc.eq.1) then                                               
         vn=alpha(n)                                                       
      else
         vn=beta(n)
      end if
c
c     *** apparent velocity ***
c
      c=vn/sin(ang)                                                     
c
      do 2 m=1,n                                                        
c
        cal=real((c/alpha(m))**2-1. )                                   
        if (cal.ge.0.) then
           ralpha(m)=csqrt((c/alpha(m))**2-1. )                            
        else
           ralpha(m)=-imag*csqrt(1. -(c/alpha(m))**2)
        end if                      
c
        cbe=real((c/beta(m))**2-1. )                                    
        if (cbe.ge.0.) then
           rbeta(m)=csqrt((c/beta(m))**2-1. )                              
        else
           rbeta(m)=-imag*csqrt(1. -(c/beta(m))**2)
        end if                        
c
        gamma(m)=2.*(beta(m)/c)**2                                      
c
    2 continue                                                          
c
      if (lc.le.2) then
c                                                                       
c     *** inverse matrix at the lowermost interface for P-SV ***
c
         en(1,1)=-2.*(beta(n)/alpha(n))**2                                 
         en(1,2)=0.                                                        
         en(1,3)=1./(rho(n)*alpha(n)**2)                                   
         en(1,4)=0.                                                        
         en(2,1)=0.                                                        
         en(2,2)=c**2*(gamma(n)-1.)/(alpha(n)**2*ralpha(n))                
         en(2,3)=0.                                                        
         en(2,4)=1./(rho(n)*alpha(n)**2*ralpha(n))                         
         en(3,1)=(gamma(n)-1.)/(gamma(n)*rbeta(n))                         
         en(3,2)=0.                                                        
         en(3,3)=-1./(rho(n)*c**2*gamma(n)*rbeta(n))                       
         en(3,4)=0.                                                        
         en(4,1)=0.                                                        
         en(4,2)=1.                                                        
         en(4,3)=0.                                                        
         en(4,4)=1./(rho(n)*c**2*gamma(n))                                 
         coef1=2.*c**2/alpha(n)**2                                         
         coef2=coef1/ralpha(n)                                             
         coef3=2.0/gamma(n)                                                
         coef4=coef3/rbeta(n)                                              
c
      end if
c
c     *** loop on frequency ***
c
      do 5 kk=1,nb2                                                     
c
        nn=n-1                                                          
        do 3 m=1,nn                                                     
          p(m)=w(kk)*ralpha(m)*h(m)*rs(m)/c                             
          q(m)=w(kk)*rbeta(m)*h(m)*rs(m)/c                              
          ! Intrinsic math functions are overloaded for complex input arguments,
          ! see https://gcc.gnu.org/onlinedocs/gcc-4.9.4/gfortran.pdf, Section 8.
          cp = cos(p(m))
          sp = sin(p(m))
          cq = cos(q(m))
          sq = sin(q(m))
          bm=beta(m)**2*rho(m)*rbeta(m)/rs(m)

c         **  matrix at layer interfaces  **                             
c
          if (lc.le.2) then
c
c            ** for P-SV problem
c
             a(1,1,m)=gamma(m)*cp-(gamma(m)-1.0)*cq                        
             a(1,2,m)=imag*((gamma(m)-1.0)*sp/ralpha(m)                    
     &            +gamma(m)*rbeta(m)*sq)                                
             a(1,3,m)=-(cp-cq)/(rho(m)*c**2)                               
             a(1,4,m)=imag*(sp/ralpha(m)+rbeta(m)*sq)/(rho(m)*c**2)        
             a(2,1,m)=-imag*(gamma(m)*ralpha(m)*sp                         
     &            +(gamma(m)-1.0)*sq/rbeta(m))                          
             a(2,2,m)=-(gamma(m)-1.0)*cp+gamma(m)*cq                       
             a(2,3,m)=imag*(ralpha(m)*sp+sq/rbeta(m))/(rho(m)*c**2)        
             a(2,4,m)=a(1,3,m)                                             
             a(3,1,m)=rho(m)*c**2*gamma(m)*(gamma(m)-1.0)*(cp-cq)          
             a(3,2,m)=imag*(rho(m)*c**2*((gamma(m)-1.0)**2*sp/ralpha(m)    
     &            +gamma(m)**2*rbeta(m)*sq))                            
             a(3,3,m)=a(2,2,m)                                             
             a(3,4,m)=a(1,2,m)                                             
             a(4,1,m)=imag*(rho(m)*c**2*(gamma(m)**2*ralpha(m)*sp          
     &            +(gamma(m)-1.0)**2*sq/rbeta(m)))                      
             a(4,2,m)=a(3,1,m)                                             
             a(4,3,m)=a(2,1,m)                                             
             a(4,4,m)=a(1,1,m)                                             
          else
c
c            ** for SH problem **
c
             al(1,1,m)=cq                                                  
             al(2,2,m)=al(1,1,m)                                           
             al(1,2,m)=imag*sq/bm                                          
             al(2,1,m)=imag*sq*bm
c
          end if
c                                          
    3   continue                                                        
c
        if (lc.le.2) then
c
c       **  matrix product for P-SV case **
c
           m=n-1                                                           
    4      continue                                                        
           if(m.eq.1)go to 40                                              

           do i=1,4
             do k=1,4
               s=0.
               do j=1,4
                 s=s+a(i,j,m)*a(j,k,m-1)
               enddo
               b(i,k,m-1)=s
             enddo
           enddo

           do i=1,4
             do k=1,4
               a(i,k,m-1)=b(i,k,m-1)
             enddo
           enddo

           m=m-1                                                           
           go to 4                                                         
   40      continue                                                        

           do i=1,4
             do k=1,4
               s=0.
               do j=1,4
                 s=s+en(i,j)*a(j,k,1)
               enddo
               aj(i,k)=s
             enddo
           enddo

           da=(aj(1,1)-aj(2,1))*(aj(3,2)-aj(4,2))-(aj(1,2)-aj(2,2))        
     &          *(aj(3,1)-aj(4,1))                                            
        if (lc.eq.1) then
c
c       **  up horizontal/total, wp vertical/total of P  **             
c
           dc=aj(4,2)-aj(3,2)                                              
           de=aj(4,1)-aj(3,1)                                              
           dca=dc/da                                                       
           dea=de/da                                                       
           up(kk)=dca*coef1*alpha(n)/c                                     
           wp(kk)=dea*coef2*alpha(n)*ralpha(n)/c
c                           
        else if (lc.eq.2) then
c
c     **  usv horizontal/total, wsv vertical/total of SV  **          
c
           gc=aj(1,2)-aj(2,2)                                              
           ge=aj(2,1)-aj(1,1)                                              
           gca=gc/da                                                       
           gea=ge/da                                                       
           usv(kk)=gca*coef4*beta(n)*rbeta(n)/c                            
           wsv(kk)=gea*coef3*beta(n)/c                                     
c
        end if
c
      else
c
c       **  matrix product for SH waves  **                             
c
         m=n-1                                                           
  114    continue                                                        
         if(m.eq.1) go to 140                                            

         do i=1,2
           do k=1,2
             sl=0.0
             do j=1,2
               sl=sl+al(i,j,m)*al(j,k,m-1)
             enddo
             bl(i,k,m-1)=sl
           enddo
         enddo

         do i=1,2
           do k=1,2
             al(i,k,m-1)=bl(i,k,m-1)
           enddo
         enddo

         m=m-1                                                           
         go to 114                                                       
  140    continue                                                        
         bn=beta(n)**2*rho(n)*rbeta(n)/rs(n)                             
         ab=al(2,1,1)+bn*al(1,1,1)
c                                       
c       **  vsh horizontal/total of SH  **                              
c
         rvs=2.*bn/ab                                                   
         vsh(kk)=rvs
c
      end if
c                                                     
    5 continue                                                          
c 
      return                                                            
      end                                                               
