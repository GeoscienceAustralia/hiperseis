      subroutine ydiz(clats,clons,clatr,clonr,
     &                delta,deltak,deltar)
cc----------------------------------------------------------------------
c Calculates distance
c for spheroidal earth between
c specified geographic source and
c receiver station coordinates
c YDAZ modified 2015
cc----------------------------------------------------------------------
c+ copyright B.L.N. KENNETT
c+ R.S.E.S. A.N.U. January 1978
cc---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real clats,clons,clatr,clonr,delta,deltak,deltar
Cf2py intent(in) clats,clons,clatr,clonr
Cf2py intent(out) delta, deltak, deltar
c radius on spheroid
      gra(x,y,e) = dsqrt( (1.0d0-e)**2 /
     &                  ((1.0d0-e*y)**2 + e*e*x*y ) )
      ecc = 0.003367
      re = 6378.388
      ec1 = (1.0d0-ecc)**2
      pi = 3.141592653589793
      pib2 = pi/2.0
      degr = pi/180.0
      dlats = clats*degr
      dlons = clons*degr
      dlatr = clatr*degr
      dlonr = clonr*degr
c geocentric coordinates
      aa=ec1*sin(dlats)
      bb=cos(dlats)
      glats = datan2 (aa,bb)
      glatr = datan2 ( ec1*sin(dlatr) ,1.0d0*cos(dlatr) )
      sps = sin(glats)**2
      cps = cos(glats)**2
      spr = sin(glatr)**2
      cpr = cos(glatr)**2
c radii at source,receiver
      rs = re*gra(sps,cps,ecc)
      rr = re*gra(spr,cpr,ecc)
c
      trs = pib2 - glats
      prs = dlons
      trr = pib2 - glatr
      prr = dlonr
c direction cosines for source
      AS = dsin(trs)*dcos(prs)
      BS = dsin(trs)*dsin(prs)
      CS = dcos(trs)
c direction cosines for receiver
      AR = dsin(trr)*dcos(prr)
      BR = dsin(trr)*dsin(prr)
      CR = dcos(trr)
c djstance
      cosdr = AS*AR + BS*BR + CS*CR
      deltar = dacos(cosdr)
c
      deltak = deltar*0.5d0*(rr+rs)
      delta = deltar/degr
c
      return
      end
