subroutine get_units(scale)
 use parameters
 use derived_types

 type (units) :: scale

 scale%length  = 1.496d13
 scale%time    =(3600d0*365.25d0*24d0)/(two*pi)
 scale%density =one/(gcgs*scale%time**2)
 scale%mass    =scale%density*scale%length**3
 scale%vel     =scale%length/scale%time
 scale%eps     =scale%vel**2*scale%mass/(scale%length)**3
 scale%kelvin  = one
 scale%rgas    = rgasCGS/(scale%vel**2/(scale%kelvin))

end subroutine
 
