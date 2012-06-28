module grid_commons
 use parameters
 use derived_types
 integer::nhex,ngrid,nbound,id_grid_center,ccol,nghost,nanchor

 integer,dimension(:),allocatable::indx_bound,column_number,indx_ghost,indx_anchor

 real(pre)::dt,time
 real(pre),dimension(:),allocatable::p,adindx,muc_array
 real(pre),dimension(:),allocatable::phi,rhotot
 real(pre),dimension(:,:),allocatable::u,qq,gforce,pforce
 real(pre),dimension(:,:),allocatable,target::cons_old,cons_new,cons
 real(pre),dimension(:,:),pointer::cons_pt
 real(pre),dimension(:),allocatable::tk_table,rho_table,deng_eos_array
 real(pre),dimension(:,:),allocatable::gamma_table2,eng_table2,gamma_table,eng_table,p_table
 real(pre),dimension(:,:),allocatable::muc_table,muc_table2,tk_table2

! the following variables need to be used before a parallel region
! is declared, and must be given here as a result.
 real(pre)::max_den_change_old,max_den_old,max_den,dtrate,maxv,maxerr
 real(pre)::muc,rate_expand=zero


 type (gridcell),dimension(:),allocatable::grid
end module grid_commons
 
