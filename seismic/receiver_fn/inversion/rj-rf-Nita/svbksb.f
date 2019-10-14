C Tiny Revision by Sheng due to occasionally strange memory accessing
C Sheng Wang
C Aug-2018 RSES ANU
      subroutine svbksb(u,w,v,m,n,mp,np,b,x)
      implicit none
      integer m,mp,n,np,NMAX
      parameter (NMAX=1000)
      integer idx,idx_j,idx_k
      double precision b(mp),u(mp,np),v(np,np),w(np),x(np)
      double precision s,tmp(NMAX)
      do 12 idx_j=1,n
          s=0.
          if(w(idx_j).ne.0.)then
              do 11 idx=1,m
                  s=s+u(idx,idx_j)*b(idx)
11            continue
              s=s/w(idx_j)
          endif
          tmp(idx_j)=s
12    continue
      do 14 idx_j=1,n
          s=0.
          do 13 idx_k=1,n
              s = s+v(idx_j,idx_k)*tmp(idx_k)
13        continue
          x(idx_j)=s
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
