! Re-write by Sheng due to occasionally strange memory accessing
! Sheng Wang
! Aug-2018 RSES ANU
subroutine voro2qmodel(voro,nmod,nmod_max,dep_min,dep_max,vs,thickness,vpvs,qa,qb)
implicit none
    real voro(nmod_max,3),dep_min,dep_max,vs(nmod_max)
    real vpvs(nmod_max),qa(nmod_max),qb(nmod_max),maxx,minn,tmp_sum,thickness(nmod_max)
    integer nmod,idx,idx_i,idx_j,order(nmod),nmod_max

    vs   = 0
    vpvs = 0
    qa   = 0
    qb   = 0
    thickness = 0
    idx  = 1

    do idx_i = 1,nmod
        maxx = dep_max
        if (idx_i == 1) then
            minn = dep_min
        else
            minn = voro(order(idx_i-1), 1)
        endif
        do idx_j = 1,nmod
            if ( (minn<voro(idx_j,1) ).and.( voro(idx_j,1)<maxx) ) then
                idx  = idx_j
                maxx = voro(idx_j,1)
            endif
        enddo
        order(idx_i) = idx
    enddo

    tmp_sum = 0.0
    do idx_i=1,nmod-1
        thickness(idx_i) = (voro(order(idx_i), 1)+voro(order(idx_i+1), 1))
        thickness(idx_i) = thickness(idx_i) * 0.5 - tmp_sum
        tmp_sum = tmp_sum + thickness(idx_i)
        vs(idx_i) = voro(order(idx_i),2)
        vpvs(idx_i) = voro(order(idx_i),3)

        qa(idx_i)= 1450
        qb(idx_i)= 600
    enddo

    thickness(nmod)=0
    vs(nmod) = voro(order(nmod),2)
    vpvs(nmod) = voro(order(nmod),3)
    qa(nmod)= 1450
    qb(nmod)= 600

    ! check correctness of generated model
    do idx_i = 1,nmod
        if ( (thickness(idx_i)>dep_max) .or. (thickness(idx_i)<dep_min) ) then
            write(*, *) "Err: Wrong thickness settings for model:"
            write(*, *) thickness
            stop
        endif
        if ( (vs(idx_i)>6.0) .or. (vs(idx_i)<0.0) ) then
            write(*, *) "Err: Wrong vs settings for model:"
            write(*, *) vs
            stop
        endif
        if ( (vpvs(idx_i)>2.0) .or. (vpvs(idx_i)<1.0) ) then
            write(*, *) "Err: Wrong vpvs settings for model:"
            write(*, *) vpvs
            stop
        endif
    enddo
    return
end