  if(convergence)then
     if(mpi_master)write(*,"(A,ES15.7)")bold_green("error="),err
  else
     if(err < eps)then
        if(mpi_master)write(*,"(A,ES15.7)")bold_yellow("error="),err
     else
        if(mpi_master)write(*,"(A,ES15.7)")bold_red("error="),err
     endif
  endif

