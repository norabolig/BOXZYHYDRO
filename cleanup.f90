subroutine cleanup
 use parameters
 use grid_commons
 implicit none
 
 deallocate(grid,phi,p,adindx,cons,cons_old,cons_new,u,qq)

end subroutine

