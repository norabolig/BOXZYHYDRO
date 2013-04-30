subroutine cleanup
 use parameters
 use grid_commons
 use derived_types
 implicit none
 
 if(allocated(grid))deallocate(grid)
 if(allocated(phi))deallocate(phi)
 if(allocated(p))deallocate(p)
 if(allocated(adindx))deallocate(adindx)
 if(allocated(cons))deallocate(cons)
 if(allocated(cons_old))deallocate(cons_old)
 if(allocated(u))deallocate(u)

end subroutine

