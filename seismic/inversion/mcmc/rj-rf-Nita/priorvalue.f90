subroutine priorvalue(point,d_min,d_max,beta_min,beta_max,width,pv_min,pv_max)
implicit none
real point,d_min,d_max,beta_min,beta_max,width,pv_min,pv_max,v

v=(beta_max-beta_min)*(point-d_min)/(d_max-d_min)+beta_min
pv_min=v-width
pv_max=v+width

! if (point<0) then
! pv_min = beta_min-width
! pv_max = beta_max+width
! endif


return
end