! Re-write by Sheng due to occasionally strange memory accessing
! Sheng Wang
! Aug-2018 RSES ANU
subroutine whichcell(point,voro,nmod,nmod_max,idx)
implicit none
    real    point,voro(nmod_max,3)
    integer nmod,idx,idx_i,nmod_max

    idx=1
    do idx_i=1,nmod
        if (abs(point-voro(idx_i,1))<abs(point-voro(idx,1))) then
            idx=idx_i
        endif
    enddo
    return
end