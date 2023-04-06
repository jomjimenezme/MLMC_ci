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

void rhs_define( vector_double rhs, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading here
  if(threading->thread != 0)
    return;

  int start = threading->start_index[l->depth];
  int end = threading->end_index[l->depth];

  if ( g.rhs == 0 ) {
    vector_double_define( rhs, 1, start, end, l );
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = ones\n");
    END_MASTER(threading)
  } else if ( g.rhs == 1 )  {
    vector_double_define( rhs, 0, start, end, l );
    if ( g.my_rank == 0 ) {
      START_LOCKED_MASTER(threading)
      rhs[0] = 1.0;
      END_LOCKED_MASTER(threading)
    }
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = first unit vector\n");
    END_MASTER(threading)
  } else if ( g.rhs == 2 ) {
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = random\n");
    END_MASTER(threading)
  } else if ( g.rhs == 3 ) {
    vector_double_define( rhs, 0, start, end, l );
  } else {
    ASSERT( g.rhs >= 0 && g.rhs <= 4 );
  }
    
}


int wilson_driver( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  int iter = 0, start = threading->start_index[l->depth], end = threading->end_index[l->depth];
  
  vector_double rhs = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.b:g.p.b;
  vector_double sol = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.x:g.p.x;

#ifdef WILSON_BENCHMARK
  START_MASTER(threading)
  prof_init( l );
  END_MASTER(threading)
  double t = -MPI_Wtime();
  double t_min = 1000;
  for ( int i=0; i<100; i++ ) {
    double tmp_t = -MPI_Wtime();
#endif
  
  vector_double_copy( rhs, source, start, end, l );  
  if ( g.method == -1 ) {
    cgn_double( &(g.p), l, threading );
  } else if ( g.mixed_precision == 2 ) {

    //double t0x=0, t1x=0, elap_time=0;
    //t0x = MPI_Wtime();

    iter = fgmres_MP( &(g.p_MP), l, threading );

    //t1x = MPI_Wtime();
    //elap_time = t1x-t0x;
    //if (g.my_rank==0) printf("elapsed time (solve phase): %-8.4lf seconds\n", elap_time);

  } else {
    iter = fgmres_double( &(g.p), l, threading );
  }
  vector_double_copy( solution, sol, start, end, l );
#ifdef WILSON_BENCHMARK
    tmp_t += MPI_Wtime();
    if ( tmp_t < t_min )
      t_min = tmp_t;
  }
  t +=MPI_Wtime();
  START_MASTER(threading)
  printf0("average over 100 solves: %lf seconds\n", t/100 );
  printf0("minimum out of 100 solves: %lf seconds\n", t_min );
  prof_print( l );
  END_MASTER(threading)
#endif
  
  return iter;
}


void solve( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  if ( g.vt.evaluation ) {
    vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    scan_var( &(g.vt), l );
    END_LOCKED_MASTER(threading)
  } else {
    wilson_driver( solution, source, l, threading );
  }
}


void solve_driver( level_struct *l, struct Thread *threading ) {
  
  vector_double solution = NULL, source = NULL;
  double minus_twisted_bc[4], norm;
 
  if(g.bc==2)
    for ( int i=0; i<4; i++ )
      minus_twisted_bc[i] = -1*g.twisted_bc[i];
  
#ifdef HAVE_TM1p1
  if( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) { 
    data_layout_n_flavours( 2, l, threading );
    printf0("inverting doublet operator\n");
  }
#endif
  PUBLIC_MALLOC( solution, complex_double, l->inner_vector_size );
  PUBLIC_MALLOC( source, complex_double, l->inner_vector_size );

  rhs_define( source, l, threading );

  if(g.bc==2)
    apply_twisted_bc_to_vector_double( source, source, g.twisted_bc, l);

  norm = global_norm_double( source, 0, l->inner_vector_size, l, threading );
  printf0("source vector norm: %le\n",norm);

#ifdef HAVE_TM1p1
  if( g.n_flavours == 1 )
#endif
#ifdef HAVE_TM
  if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
    if(g.downprop) {
      
      START_MASTER(threading)  
      printf0("\n\n+--------------------------- up ---------------------------+\n\n");
      END_MASTER(threading)

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

      // calling the coarsest-level solver once on setup
//#if defined(GCRODR) || defined(POLYPREC) || defined(BLOCK_JACOBI)
#if defined(GCRODR)
      {
        level_struct *lx = l;

        START_MASTER(threading)
        printf0( "\nPre-constructing coarsest-level data ...\n" );
        END_MASTER(threading)

        while (1) {
          if ( lx->level==0 ) {

            if ( !(lx->idle) ) {

            if ( g.mixed_precision==0 ) {

              gmres_double_struct* px = &(lx->p_double);

              // set RHS to random
              START_MASTER(threading)
              vector_double_define_random( px->b, px->v_start, px->v_end, lx );
              END_MASTER(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 1;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              double buff1x = px->tol;
              double buff2x = g.coarse_tol;
              START_MASTER(threading)
              px->tol = 1.0e-20;
              g.coarse_tol = 1.0e-20;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)
              // call the coarsest-level solver
              while ( px->gcrodr_double.CU_usable==0 ) {
                coarse_solve_odd_even_double( px, &(lx->oe_op_double), lx, threading );
              }
              START_MASTER(threading)
              px->tol = buff1x;
              g.coarse_tol = buff2x;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 0;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

            }
            else {

              gmres_float_struct* px = &(lx->p_float);

              // set RHS to random
              START_MASTER(threading)
              vector_float_define_random( px->b, px->v_start, px->v_end, lx );
              END_MASTER(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 1;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              double buff1x = px->tol;
              double buff2x = g.coarse_tol;
              START_MASTER(threading)
              px->tol = 1.0e-20;
              g.coarse_tol = 1.0e-20;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)
              // call the coarsest-level solver
              while ( px->gcrodr_float.CU_usable==0 ) {
                coarse_solve_odd_even_float( px, &(lx->oe_op_float), lx, threading );
              }
              START_MASTER(threading)
              px->tol = buff1x;
              g.coarse_tol = buff2x;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 0;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

            }

            } // end of !idle if

            break;
          }
          else { lx = lx->next_level; }
        }

        START_MASTER(threading)
        printf0( "... done\n\n" );
        END_MASTER(threading)

      }
#endif

      START_MASTER(threading)
      g.avg_b1 = 0.0;
      g.avg_b2 = 0.0;
      g.avg_crst = 0.0;
      END_MASTER(threading)
      solve( solution, source, l, threading );
      START_MASTER(threading)
      printf0( "avg coarsest iters = %f\n",g.avg_crst );
      END_MASTER(threading)

      if(g.bc==2)
      apply_twisted_bc_to_vector_double( solution, solution, minus_twisted_bc, l);
      
      START_LOCKED_MASTER(threading)  
      printf0("\n\n+-------------------------- down --------------------------+\n\n");
      g.mu*=-1;
      g.mu_odd_shift*=-1;
      g.mu_even_shift*=-1;
      END_LOCKED_MASTER(threading)
  
      tm_term_update( g.mu, l, threading );
      finalize_operator_update( l, threading );
    } 
#endif

#ifdef POLYPREC
  {

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_MASTER(threading)

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

    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

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

      // calling the coarsest-level solver once on setup
//#if defined(GCRODR) || defined(POLYPREC) || defined(BLOCK_JACOBI)
#if defined(GCRODR)
      {
        level_struct *lx = l;

        START_MASTER(threading)
        printf0( "\nPre-constructing coarsest-level data ...\n" );
        END_MASTER(threading)

        while (1) {
          if ( lx->level==0 ) {

            if ( !(lx->idle) ) {

            if ( g.mixed_precision==0 ) {

              gmres_double_struct* px = &(lx->p_double);

              // set RHS to random
              START_MASTER(threading)
              vector_double_define_random( px->b, px->v_start, px->v_end, lx );
              END_MASTER(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 1;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              double buff1x = px->tol;
              double buff2x = g.coarse_tol;
              START_MASTER(threading)
              px->tol = 1.0e-20;
              g.coarse_tol = 1.0e-20;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)
              // call the coarsest-level solver
              while ( px->gcrodr_double.CU_usable==0 ) {
                coarse_solve_odd_even_double( px, &(lx->oe_op_double), lx, threading );
              }
              START_MASTER(threading)
              px->tol = buff1x;
              g.coarse_tol = buff2x;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 0;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

            }
            else {

              gmres_float_struct* px = &(lx->p_float);

              // set RHS to random
              START_MASTER(threading)
              vector_float_define_random( px->b, px->v_start, px->v_end, lx );
              END_MASTER(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 1;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              double buff1x = px->tol;
              double buff2x = g.coarse_tol;
              START_MASTER(threading)
              px->tol = 1.0e-20;
              g.coarse_tol = 1.0e-20;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)
              // call the coarsest-level solver
              while ( px->gcrodr_float.CU_usable==0 ) {
                coarse_solve_odd_even_float( px, &(lx->oe_op_float), lx, threading );
              }
              START_MASTER(threading)
              px->tol = buff1x;
              g.coarse_tol = buff2x;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

              START_MASTER(threading)
              g.gcrodr_calling_from_setup = 0;
              END_MASTER(threading)
              SYNC_MASTER_TO_ALL(threading)

            }

            } // end of !idle if

            break;
          }
          else { lx = lx->next_level; }
        }

        START_MASTER(threading)
        printf0( "... done\n\n" );
        END_MASTER(threading)

      }
#endif

  START_MASTER(threading)
  g.avg_b1 = 0.0;
  g.avg_b2 = 0.0;
  g.avg_crst = 0.0;
  END_MASTER(threading)

  solve( solution, source, l, threading );

  START_MASTER(threading)
  printf0( "avg coarsest iters = %f\n",g.avg_crst );
  END_MASTER(threading)

  if(g.bc==2)
    apply_twisted_bc_to_vector_double( solution, solution, minus_twisted_bc, l);

  norm = global_norm_double( solution, 0, l->inner_vector_size, l, threading );
  printf0("solution vector norm: %le\n",norm);

  PUBLIC_FREE( solution, complex_double, l->inner_vector_size );
  PUBLIC_FREE( source, complex_double, l->inner_vector_size );

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) 
    data_layout_n_flavours( 1, l, threading );
#endif
}

