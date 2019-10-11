      FUNCTION dpythag(a,b)
      DOUBLE PRECISION a,b,dpythag
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        dpythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          dpythag=0.0d0
        else
          dpythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
