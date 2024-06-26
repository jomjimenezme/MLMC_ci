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


// aux functions
void gcrodr_and_polyprec_some_resets1_float( level_struct *l );
void gcrodr_and_polyprec_some_resets1_double( level_struct *l );
void gcrodr_and_polyprec_some_resets2( level_struct *l );
void gcrodr_and_polyprec_some_resets3( level_struct *l, struct Thread *threading );


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
    g.if_rademacher=0;
    g.on_solve = 0;
    struct Thread threading;
    struct Thread *threadingx = &threading;  
    setup_threading(&threading, commonthreaddata, &l);
    setup_no_threading(no_threading, &l);

    // initial phase of AMG
    double t_setup0, t_update0, t_setup1, t_update1;
    t_setup0 = MPI_Wtime();
    method_setup( NULL, &l, &threading );
    t_setup1 = MPI_Wtime();
    START_MASTER(threadingx)
    if(g.my_rank==0)printf("Time for initial phase of DDalphaAMG : %f\n", t_setup1-t_setup0);
    END_MASTER(threadingx)
    fflush(0);

    gcrodr_and_polyprec_some_resets1_float( &l );

     // iterative phase of AMG
    t_update0 = MPI_Wtime();
    method_update( l.setup_iter, &l, &threading );
    t_update1 = MPI_Wtime();

    START_MASTER(threadingx)
    if(g.my_rank==0)printf("Time for iterative phase of DDalphaAMG : %f\n", t_update1-t_update0);
    END_MASTER(threadingx)
    fflush(0);

    gcrodr_and_polyprec_some_resets1_double( &l );
    gcrodr_and_polyprec_some_resets2( &l );
    gcrodr_and_polyprec_some_resets3( &l,threadingx );

    // calling block power iteration, to generate all deflation vectors at all levels
    g.on_solve = 1;
    double t_powerit0, t_powerit1;
    t_powerit0 = MPI_Wtime();
    block_powerit_driver_double( &l, &threading );
    t_powerit1 =MPI_Wtime();
    START_MASTER(threadingx)
    if(g.my_rank==0)printf("Time for block power iteration : %f\n", t_powerit1-t_powerit0);
    END_MASTER(threadingx)
    fflush(0);

    // pre-setting some values for the Hutchinson struct
    {
      //l.h_double.max_iters = 1000;
      //l.h_double.min_iters = 1000;

      l.h_double.trace_tol = 1.0e-4;
      hutchinson_diver_double_init( &l, &threading );  
      hutchinson_diver_double_alloc( &l, &threading );

      hutchinson_double_struct* h = &(l.h_double);
      h->tol_per_level = malloc(sizeof(double)*l.depth);
      h->tol_per_level[0] = sqrt(0.75);
      h->tol_per_level[1] = sqrt(0.20);
      h->tol_per_level[2] = sqrt(0.05);
    }

    // legacy line, from when this was a solver library only
    //solve_driver( &l, &threading );

    /*
      now, we call the method to compute the trace. This is first based on <d0 trace deflation type>
      from the input file : if this has a value less than 3, then we take this as the method of choice
      to compute the trace. If, in turn, its value is 3, then we check the value set via <d0 trace op type>
      to set the trace-computation method
    */

    complex_double trace;

    int op_type;
    if ( g.trace_deflation_type[0]==1 || g.trace_deflation_type[0]==2 ) {
      if ( g.trace_deflation_type[0] != g.trace_op_type )
        if(g.my_rank==0) warning("CAREFUL : you are using different operators for BPI and the trace method");
    } else if ( g.trace_deflation_type[0]>=3 && g.trace_deflation_type[0]<=5 ) {
      if ( g.trace_op_type<3 || g.trace_op_type>5 ) {
        if(g.my_rank==0) warning("CAREFUL : you are using different operators for BPI and the trace method");
      }
    }
    op_type = g.trace_op_type;

    if(g.my_rank==0) printf("\n");
    if ( op_type==2 ) {
      START_MASTER(threadingx)
      if(g.my_rank==0) printf("Using deflated plain Hutchinson for computing the trace\n");
      END_MASTER(threadingx)

      trace = hutchinson_driver_double( &l, &threading );

      START_MASTER(threadingx)
      if(g.my_rank==0) printf("\n");
      if(g.my_rank==0) printf("Resulting trace from deflated plain Hutchinson = %f+i%f\n", CSPLIT(trace));
      END_MASTER(threadingx)
    } else if ( op_type==1 ) {
      START_MASTER(threadingx)
      if(g.my_rank==0) printf("Using (traditional) MGMLMC for computing the trace\n");
      END_MASTER(threadingx)

      trace = mlmc_hutchinson_driver_double( &l, &threading );

      START_MASTER(threadingx)
      if(g.my_rank==0) printf("\n");
      if(g.my_rank==0) printf("Resulting trace from (traditional) MGMLMC = %f+i%f\n", CSPLIT(trace));
      END_MASTER(threadingx)
    } else if ( op_type>=3 && op_type<=5 ) {
      START_MASTER(threadingx)
      if(g.my_rank==0) printf("Using (split) MGMLMC for computing the trace\n");
      END_MASTER(threadingx)

      trace = split_mlmc_hutchinson_driver_double( &l, &threading );

      START_MASTER(threadingx)
      if(g.my_rank==0) printf("\n");
      if(g.my_rank==0) printf("Resulting trace from (split) MGMLMC = %f+i%f\n", CSPLIT(trace));
      END_MASTER(threadingx)
    }

    hutchinson_diver_double_free( &l, &threading );
    // FIXME : uncommenting this gives seg fault
    //block_powerit_double_free( &l, &threading );
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









// -------------------------------------------------------------------------------------


// aux functions


void gcrodr_and_polyprec_some_resets1_float( level_struct *l ){
#if defined(POLYPREC) || defined(GCRODR)
  {
    level_struct *lx = l;
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
}


void gcrodr_and_polyprec_some_resets1_double( level_struct *l ){
#if defined(POLYPREC) || defined(GCRODR)
  {
    level_struct *lx = l;
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
}


void gcrodr_and_polyprec_some_resets2( level_struct *l ){
#ifdef POLYPREC
  {
    // setting flag to re-update lejas
    level_struct *lx = l;
    while (1) {
      if ( lx->level==0 ) {
        if ( g.mixed_precision==0 ) {
          lx->p_double.polyprec_double.update_lejas = 1;
          lx->p_double.polyprec_double.preconditioner = NULL;
        }
        else {
          lx->p_float.polyprec_float.update_lejas = 1;
          lx->p_float.polyprec_float.preconditioner = NULL;
        }
        break;
      }
      else { lx = lx->next_level; }
    }
  }
#endif

#ifdef BLOCK_JACOBI
  {
    // setting flag to re-update lejas
    level_struct *lx = l;
    while (1) {
      if ( lx->level==0 ) {
        if ( g.mixed_precision==0 ) {
          lx->p_double.block_jacobi_double.local_p.polyprec_double.update_lejas = 1;
          lx->p_double.block_jacobi_double.BJ_usable = 0;
        }
        else {
          lx->p_float.block_jacobi_float.local_p.polyprec_float.update_lejas = 1;
          lx->p_float.block_jacobi_float.BJ_usable = 0;
        }
        break;
      }
      else { lx = lx->next_level; }
    }
  }
#endif

#ifdef GCRODR
  {
    // setting flag to re-update recycling subspace
    level_struct *lx = l;
    while (1) {
      if ( lx->level==0 ) {
        if ( g.mixed_precision==0 ) {
          //lx->p_double.gcrodr_double.CU_usable = 0;
          lx->p_double.gcrodr_double.update_CU = 1;
          lx->p_double.gcrodr_double.upd_ctr = 0;
        }
        else {
          //lx->p_float.gcrodr_float.CU_usable = 0;
          lx->p_float.gcrodr_float.update_CU = 1;
          lx->p_float.gcrodr_float.upd_ctr = 0;
        }
        break;
      }
      else { lx = lx->next_level; }
    }
  }
#endif
}


void gcrodr_and_polyprec_some_resets3( level_struct *l, struct Thread *threadingx ){
  // calling the coarsest-level solver once on setup
#if defined(GCRODR)
  {
    level_struct *lx = l;

    START_MASTER(threadingx)
    printf0( "\nPre-constructing coarsest-level data ...\n" );
    END_MASTER(threadingx)

    while (1) {
      if ( lx->level==0 ) {

        if ( !(lx->idle) ) {

        if ( g.mixed_precision==0 ) {

          gmres_double_struct* px = &(lx->p_double);

          // set RHS to random
          START_MASTER(threadingx)
          vector_double_define_random( px->b, px->v_start, px->v_end, lx );
          END_MASTER(threadingx)

          START_MASTER(threadingx)
          g.gcrodr_calling_from_setup = 1;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)

          double buff1x = px->tol;
          double buff2x = g.coarse_tol;
          START_MASTER(threadingx)
          px->tol = 1.0e-20;
          g.coarse_tol = 1.0e-20;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)
          // call the coarsest-level solver
          while ( px->gcrodr_double.CU_usable==0 ) {
            coarse_solve_odd_even_double( px, &(lx->oe_op_double), lx, threadingx );
          }
          START_MASTER(threadingx)
          px->tol = buff1x;
          g.coarse_tol = buff2x;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)

          START_MASTER(threadingx)
          g.gcrodr_calling_from_setup = 0;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)

        }
        else {

          gmres_float_struct* px = &(lx->p_float);

          // set RHS to random
          START_MASTER(threadingx)
          vector_float_define_random( px->b, px->v_start, px->v_end, lx );
          END_MASTER(threadingx)

          START_MASTER(threadingx)
          g.gcrodr_calling_from_setup = 1;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)

          double buff1x = px->tol;
          double buff2x = g.coarse_tol;
          START_MASTER(threadingx)
          px->tol = 1.0e-20;
          g.coarse_tol = 1.0e-20;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)
          // call the coarsest-level solver
          while ( px->gcrodr_float.CU_usable==0 ) {
            coarse_solve_odd_even_float( px, &(lx->oe_op_float), lx, threadingx );
          }
          START_MASTER(threadingx)
          px->tol = buff1x;
          g.coarse_tol = buff2x;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)

          START_MASTER(threadingx)
          g.gcrodr_calling_from_setup = 0;
          END_MASTER(threadingx)
          SYNC_MASTER_TO_ALL(threadingx)

        }

        } // end of !idle if

        break;
      }
      else { lx = lx->next_level; }
    }

    START_MASTER(threadingx)
    printf0( "... done\n\n" );
    END_MASTER(threadingx)

  }
#endif
}
