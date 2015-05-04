!
! Define variables for global memory
!
module grid_commons
 use parameters
 use derived_types
 integer::nhex,ngrid,nbound,id_grid_center,ccol,nghost,nanchor

 integer,dimension(:),allocatable::indx_bound,column_number,indx_ghost,indx_anchor

 real(pre)::dt,time
 real(pre),dimension(:),allocatable::p,adindx,muc_array
 real(pre),dimension(:),allocatable::phi,phi_new,rhotot
 real(pre),dimension(:,:),allocatable::u,qq,gforce,pforce
 real(pre),dimension(:,:),allocatable::cons_old,cons,fluxtmp
 real(pre),dimension(:,:),allocatable::state_p_p,state_p_m
 real(pre),dimension(:,:),allocatable::state_d_p,state_d_m
 real(pre),dimension(:,:),allocatable::state_e_p,state_e_m
 real(pre),dimension(:,:),allocatable::state_g_p,state_g_m
 real(pre),dimension(:,:),allocatable::state_rmom_p,state_rmom_m
 real(pre),dimension(:,:),allocatable::state_amom_p,state_amom_m
 real(pre),dimension(:,:),allocatable::slope_p,slope_d,slope_e,slope_g
 real(pre),dimension(:,:),allocatable::slope_rmom,slope_amom
 real(pre),dimension(:,:,:),allocatable::state_u_p,state_u_m,slope_u,f_cor

 real(pre),dimension(:),allocatable::tk_table,rho_table,deng_eos_array
 real(pre),dimension(:,:),allocatable::gamma_table2,eng_table2,gamma_table,eng_table,p_table
 real(pre),dimension(:,:),allocatable::muc_table,muc_table2,tk_table2

! the following variables need to be used before a parallel region
! is declared, and must be given here as a result.
 real(pre)::max_den_change_old,max_den_old,max_den,dtrate,maxv,maxerr
 real(pre)::muc,rate_expand=zero


 type (gridcell),dimension(:),allocatable::grid
end module grid_commons
 
