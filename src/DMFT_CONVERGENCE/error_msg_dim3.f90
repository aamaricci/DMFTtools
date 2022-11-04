 if(convergence)then
     if(mpi_master)then
        write(*,"(A,ES15.7,I8)")bold_green("max error="),error(1),mpi_rank
        write(*,"(A,ES15.7,I8)")bold_green("    error="),err,mpi_rank
        write(*,"(A,ES15.7,I8)")bold_green("min error="),error(2),mpi_rank
     endif
  else
     if(err < eps)then
        if(mpi_master)then
           write(*,"(A,ES15.7,I8)")bold_yellow("max error="),error(1),mpi_rank
           write(*,"(A,ES15.7,I8)")bold_yellow("    error="),err,mpi_rank
           write(*,"(A,ES15.7,I8)")bold_yellow("min error="),error(2),mpi_rank
        endif
     else
        if(mpi_master)then
           write(*,"(A,ES15.7,I8)")bold_red("max error="),error(1),mpi_rank
           write(*,"(A,ES15.7,I8)")bold_red("    error="),err,mpi_rank
           write(*,"(A,ES15.7,I8)")bold_red("min error="),error(2),mpi_rank
        endif
     endif
  endif
