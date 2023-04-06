/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */
 
#include "main.h"

extern global_struct g;
#ifdef HAVE_HDF5
Hdf5_fileinfo h5info;
#endif
extern struct common_thread_data *commonthreaddata;
extern struct Thread *no_threading;

int main( int argc, char **argv ) {
    
#ifdef HAVE_HDF5
  h5info.filename=NULL;
  h5info.file_id=-1; 
  h5info.rootgroup_id=-1; 
  h5info.configgroup_id=-1;
  h5info.eigenmodegroup_id=-1;
  h5info.thiseigenmodegroup_id=-1;
  h5info.isOpen=0;
  h5info.mode=-1;
#endif
  level_struct l;
  config_double hopp = NULL;
  
  MPI_Init( &argc, &argv );
  
  predefine_rank( MPI_COMM_WORLD );
  if ( g.my_rank == 0 ) {
    printf("\n\n+----------------------------------------------------------+\n");
    printf("| The DDalphaAMG solver library.                           |\n");
    printf("| Copyright (C) 2016, Matthias Rottmann, Artur Strebel,    |\n");
    printf("|       Simon Heybrock, Simone Bacchio, Bjoern Leder.      |\n");
    printf("|                                                          |\n");
    printf("| This program comes with ABSOLUTELY NO WARRANTY.          |\n");
    printf("+----------------------------------------------------------+\n\n");
  }
  
  method_init( &argc, &argv, &l );
  
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);
  
  MALLOC( hopp, complex_double, 3*l.inner_vector_size );

  if(g.in_format == _LIME)
    lime_read_conf( (double*)(hopp), g.in, &(g.plaq_hopp) );
  else 
    read_conf( (double*)(hopp), g.in, &(g.plaq_hopp), &l );

  // store configuration, compute clover term
  dirac_setup( hopp, &l );
  FREE( hopp, complex_double, 3*l.inner_vector_size );

  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);
  
  THREADED(g.num_openmp_processes)
  {
    g.on_solve = 0;
    struct Thread threading;
    setup_threading(&threading, commonthreaddata, &l);
    setup_no_threading(no_threading, &l);

    //double t0x=0, t1x=0, elap_time=0;

    //t0x = MPI_Wtime();

    // setup up initial MG hierarchy
    method_setup( NULL, &l, &threading );

    //t1x = MPI_Wtime();
    //elap_time = t1x-t0x;
    //if (g.my_rank==0) printf("elapsed time (init setup phase): %-8.4lf seconds\n", elap_time);

    //t0x = MPI_Wtime();

#if defined(POLYPREC) || defined(GCRODR)
    {
      level_struct *lx = &l;
      while (1) {
        if ( lx->level==0 ) {
          if ( g.mixed_precision==0 ) {
#ifdef GCRODR
            lx->p_double.gcrodr_double.k = g.gcrodr_k_setup;
#endif
#ifdef POLYPREC
            lx->p_float.polyprec_float.d_poly = g.polyprec_d_setup;
#endif
          }
          else {
#ifdef GCRODR
            lx->p_float.gcrodr_float.k = g.gcrodr_k_setup;
#endif
#ifdef POLYPREC
            lx->p_float.polyprec_float.d_poly = g.polyprec_d_setup;
#endif
          }
          break;
        }
        else { lx = lx->next_level; }
      }
    }
#endif

    // iterative phase
    method_update( l.setup_iter, &l, &threading );

    //t1x = MPI_Wtime();
    //elap_time = t1x-t0x;
    //if (g.my_rank==0) printf("elapsed time (iterative setup phase): %-8.4lf seconds\n", elap_time);

#if defined(POLYPREC) || defined(GCRODR)
    {
      level_struct *lx = &l;
      while (1) {
        if ( lx->level==0 ) {
          if ( g.mixed_precision==0 ) {
#ifdef POLYPREC
            lx->p_double.polyprec_double.d_poly = g.polyprec_d_solve;
#endif
#ifdef GCRODR
            lx->p_double.gcrodr_double.k = g.gcrodr_k_solve;
#endif
          }
          else {
#ifdef POLYPREC
            lx->p_float.polyprec_float.d_poly = g.polyprec_d_solve;
#endif
#ifdef GCRODR
            lx->p_float.gcrodr_float.k = g.gcrodr_k_solve;
#endif
          }
          break;
        }
        else { lx = lx->next_level; }
      }
    }
#endif

    g.on_solve = 1;
    solve_driver( &l, &threading );
  }
  
  finalize_common_thread_data(commonthreaddata);
  finalize_no_threading(no_threading);
  free(commonthreaddata);
  free(no_threading);

  method_free( &l );
  method_finalize( &l );
  
  MPI_Finalize();
  
  return 0;
}
