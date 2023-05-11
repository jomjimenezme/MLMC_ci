
#include "main.h"





void block_powerit_PRECISION_init_and_alloc( int spec_type, int op_id, int depth_bp_op, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading ){


  int i;

  // access l at the right level
  level_struct* lx = l;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = lx->next_level;
  }

  lx->powerit.nr_vecs = nr_vecs;
  lx->powerit.bp_tol = bp_tol;
  lx->powerit.nr_cycles = nr_bpi_cycles;

  lx->powerit.spec_type = spec_type;
  lx->powerit.gs_buffer = NULL;
  PUBLIC_MALLOC( lx->powerit.gs_buffer, complex_PRECISION, 2*lx->powerit.nr_vecs );

  lx->powerit.vecs = NULL;
  PUBLIC_MALLOC( lx->powerit.vecs, complex_PRECISION*, lx->powerit.nr_vecs );
  lx->powerit.vecs[0] = NULL;
  PUBLIC_MALLOC( lx->powerit.vecs[0], complex_PRECISION, lx->powerit.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    lx->powerit.vecs[i] = lx->powerit.vecs[0] + i*lx->vector_size;
  }

  lx->powerit.vecs_buff1 = NULL;
  PUBLIC_MALLOC( lx->powerit.vecs_buff1, complex_PRECISION, lx->vector_size );

  lx->powerit.vecs_buff2 = NULL;
  PUBLIC_MALLOC( lx->powerit.vecs_buff2, complex_PRECISION, lx->vector_size );

  lx->powerit.vecs_buff3 = NULL;
  PUBLIC_MALLOC( lx->powerit.vecs_buff3, complex_PRECISION, lx->vector_size );

}


void block_powerit_PRECISION_free( level_struct* l, struct Thread* threading ){


  int i,j;

  for( j=0;j<g.num_levels;j++ ){
    // in case no deflation is requested
    if( g.trace_deflation_type[j]==3 ){ continue; }

    // access l at the right level
    level_struct* lx = l;
    for( i=0;i<g.num_levels;i++ ){
      if( i==j ){ break; }
      lx = lx->next_level;
    }

    PUBLIC_FREE( lx->powerit.vecs[0], complex_PRECISION, lx->powerit.nr_vecs*lx->vector_size );
    PUBLIC_FREE( lx->powerit.vecs, complex_PRECISION*, lx->powerit.nr_vecs );
    PUBLIC_FREE( lx->powerit.vecs_buff1, complex_PRECISION, lx->vector_size );
    PUBLIC_FREE( lx->powerit.vecs_buff2, complex_PRECISION, lx->vector_size );
    PUBLIC_FREE( lx->powerit.vecs_buff3, complex_PRECISION, lx->vector_size );

    PUBLIC_FREE( lx->powerit.gs_buffer, complex_PRECISION, 2*lx->powerit.nr_vecs );
  }

}


void block_powerit_driver_PRECISION( level_struct* l, struct Thread* threading ){

  int i,op_id,spec_type;

  // specify the following in the .ini input file, at different levels
  // dx trace deflation type: 0   // 0 is difference, 1 is non-difference, 2 is split orthogonal, 3 is no deflation
  // dx trace deflation nr vectors: 10

  for( i=0;i<g.num_levels;i++ ){
              if(g.my_rank==0)printf("\nLevel and request--------------- %d request: %d\n", i, g.trace_deflation_type[i]);

    // in case no deflation is requested
    if( g.trace_deflation_type[i]==3 ){ continue; }
    switch(g.trace_deflation_type[i]){
      case 0:
        apply_to_one_vector = powerit_diff_op_PRECISION; op_id = _NON_DIFF_OP;
        break;
      case 1:
        apply_to_one_vector = powerit_non_diff_op_PRECISION; op_id = _NON_DIFF_OP;
        break;
      case 2:
        apply_to_one_vector = powerit_split_op_PRECISION; op_id = _SPLIT_OP;
        break;

      default:
        error0("Uknown type for operator in block power iteration\n");
    }

    int depth_bp_op = i;
    int nr_bp_vecs = g.trace_deflation_nr_vectors[i];
    double bp_tol = g.trace_powerit_solver_tol[i];
    int nr_bpi_cycles = g.trace_powerit_cycles[i];

    switch(g.trace_powerit_spectrum_type[i]){
      case 0:
        spec_type = _EVs;
        break;
      case 1:
        spec_type = _SVs;
        break;
      default:
        error0("Uknown type of spectrum to be extracted\n");
    }

    // IMPORTANT :
    //		   -- always call this operation with the finest-level l
    //		   -- after calling power iteration, the result is in lx->powerit.vecs, with lx
    //		      the level struct of the chosen level
    if( depth_bp_op==(g.num_levels-1) && apply_to_one_vector == powerit_diff_op_PRECISION){
      error0("There is no difference level operator at the coarsest level\n");
    }

    block_powerit_PRECISION_init_and_alloc( spec_type, op_id, depth_bp_op, nr_bp_vecs, nr_bpi_cycles, bp_tol, l, threading );
    //blind_bp_op_PRECISION_apply( depth_bp_op, l, threading );
    block_powerit_PRECISION( op_id, depth_bp_op, l, threading );
    printf("About to finish powerit_driver\n");  
  }
printf("EXITING powerit driver\n");
}



void block_powerit_PRECISION( int op_id, int depth_bp_op, level_struct *l, struct Thread *threading ){

  // access l at the right level
  level_struct* lx = l;
  int i;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = lx->next_level;
  }

  // set the power iteration vectors to random
  START_LOCKED_MASTER(threading)
  vector_PRECISION_define_random( lx->powerit.vecs[0], 0, lx->powerit.nr_vecs*lx->vector_size, lx );
  END_LOCKED_MASTER(threading)
  SYNC_CORES(threading)

  for( i=0;i<lx->powerit.nr_cycles;i++ ){
    // apply the operator on the vectors ...
    blind_bp_op_PRECISION_apply( op_id, lx, threading );
    // ... and the resulting vectors are in lx->powerit.vecs
    //orthogonalize
    bp_qr_PRECISION( lx, threading );
  }

  // in the SVs case, this tests the eigenvectors coming out of the Hermitian problem
  test_powerit_quality_PRECISION( op_id, lx, threading );
printf("AFTER test call, applying gamma\n");
  // apply gamma5 to the final result, if singular vectors are wanted
  if( lx->powerit.spec_type ==_SVs ){
    for( i=0;i<lx->powerit.nr_vecs;i++ ){
      if( lx->depth==0 ){
        gamma5_PRECISION( lx->powerit.vecs[i], lx->powerit.vecs[i], lx, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, lx->inner_vector_size, &startg5, &endg5, lx, threading, lx->num_lattice_site_var );
        coarse_gamma5_PRECISION( lx->powerit.vecs[i], lx->powerit.vecs[i], startg5, endg5, lx );
      }
    }
    SYNC_CORES(threading)
  }
printf("AFTER Applying Gamma call\n");
}



void bp_qr_PRECISION( level_struct* lx, struct Thread* threading ){

  gram_schmidt_PRECISION( lx->powerit.vecs, lx->powerit.gs_buffer, 0, lx->powerit.nr_vecs, lx, threading );
  gram_schmidt_PRECISION( lx->powerit.vecs, lx->powerit.gs_buffer, 0, lx->powerit.nr_vecs, lx, threading );
}


void test_powerit_quality_PRECISION( int op_id, level_struct* lx, struct Thread* threading ){

  int i, start, end;
  vector_PRECISION* vecs_buff1;
  vector_PRECISION* vecs_buff2;

  gmres_PRECISION_struct* px;
  if( lx->depth==0 ){ px = &(g.p); } else { px = &(lx->p_PRECISION); }
  compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

  vecs_buff1 = NULL;
  vecs_buff2 = NULL;
  
  PUBLIC_MALLOC( vecs_buff2, complex_PRECISION*, lx->powerit.nr_vecs );
  PUBLIC_MALLOC( vecs_buff1, complex_PRECISION*, lx->powerit.nr_vecs );
  
  vecs_buff1[0] = NULL;
  vecs_buff2[0] = NULL;
  
  PUBLIC_MALLOC( vecs_buff1[0], complex_PRECISION, lx->powerit.nr_vecs*lx->vector_size );
  PUBLIC_MALLOC( vecs_buff2[0], complex_PRECISION, lx->powerit.nr_vecs*lx->vector_size );
  
  START_MASTER(threading)
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    vecs_buff1[i] = vecs_buff1[0] + i*lx->vector_size;
    vecs_buff2[i] = vecs_buff2[0] + i*lx->vector_size;
  }
  END_MASTER(threading)
  SYNC_CORES(threading)


  // backup of lx->powerit.vecs
  for( i=0;i<lx->powerit.nr_vecs;i++ ){
    vector_PRECISION_copy( vecs_buff1[i], lx->powerit.vecs[i], start, end, lx );
  }

  // apply the operator
  blind_bp_op_PRECISION_apply( op_id, lx, threading );

  // swap pointers  
  complex_PRECISION** buff_ptr = lx->powerit.vecs;
  complex_PRECISION** buff_ptr1 = vecs_buff1;    
  START_MASTER(threading)
    lx->powerit.vecs = buff_ptr1;
    vecs_buff1 = buff_ptr;
  END_MASTER(threading)  
  SYNC_MASTER_TO_ALL(threading)
  
  for( i=0;i<lx->powerit.nr_vecs;i++ ){
    // compute the Rayleigh quotient
    complex_PRECISION rq;
    rq = global_inner_product_PRECISION( lx->powerit.vecs[i], vecs_buff1[i], px->v_start, px->v_end, lx, threading );
    double norm = global_norm_PRECISION( lx->powerit.vecs[i], 0, lx->inner_vector_size, lx, threading );
    rq /= norm;

    // print the Rayleigh quotient
    START_MASTER(threading)
    if (g.my_rank==0) printf( "Rayleigh quotient = %.16f+i%.16f \t",CSPLIT(rq) );
    END_MASTER(threading)

    // compute the eigenvalue residual
    vector_PRECISION_scale( vecs_buff2[i], lx->powerit.vecs[i], rq, start, end, lx );
    vector_PRECISION_minus( vecs_buff1[i], vecs_buff1[i], vecs_buff2[i], start, end, lx );
    double resx = global_norm_PRECISION( vecs_buff1[i], 0, lx->inner_vector_size, lx, threading ) / ( cabs(rq)*norm );
    
    // print the residuals
    START_MASTER(threading)
    if (g.my_rank==0) printf( "Eigenvalue residual = %.16f\n",resx );
    END_MASTER(threading)
  }
  
  //Restore backup
  for( i=0;i<lx->powerit.nr_vecs;i++ ){
    vector_PRECISION_copy( lx->powerit.vecs[i], vecs_buff1[i], start, end, lx );
  }
  //swap pointers Back
  START_MASTER(threading)
    vecs_buff1 = buff_ptr1;
    lx->powerit.vecs = buff_ptr;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  PUBLIC_FREE( vecs_buff1[0], complex_PRECISION, lx->powerit.nr_vecs*lx->vector_size );
  PUBLIC_FREE( vecs_buff2[0], complex_PRECISION, lx->powerit.nr_vecs*lx->vector_size );
}





void powerit_non_diff_op_PRECISION( level_struct *l, int i, struct Thread *threading ){
      
  int start, end;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );

  fgmres_PRECISION( p, l, threading );
      
  vector_PRECISION_copy( l->powerit.vecs[i], p->x, start, end, l );

}



void powerit_diff_op_PRECISION( level_struct *l, int i, struct Thread *threading ){
  
  int start, end;    
  
  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );
  
      
  // coarse
  level_struct* lxc = l->next_level;
  gmres_PRECISION_struct* pxc = get_p_struct_PRECISION(lxc);

  apply_R_PRECISION(pxc->b, p->b, l, threading);     
 
  lxc->powerit.tol_buff = l->powerit.tol_buff;
  apply_solver_powerit_PRECISION(lxc, threading);

  apply_P_PRECISION(l->powerit.vecs_buff2, pxc->x, l, threading);
  
  // fine
  
  START_MASTER(threading)
  p->x = l->powerit.vecs_buff1;
  END_MASTER(threading)
  SYNC_CORES(threading)
  fgmres_PRECISION( p, l, threading );
  
  vector_PRECISION_minus( l->powerit.vecs[i], p->x, l->powerit.vecs_buff2, start, end, l );
   
}

void powerit_split_op_PRECISION( level_struct *l, int i, struct Thread *threading ){
//TODO: IMLEMENT BOTH SPLIT OPERATORS
  int start, end;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );
  //  (I - P_{l} P_{l}^{H}) A_{l}^{-1} 
  apply_solver_powerit_PRECISION(l, threading);
      
  apply_R_PRECISION(l->powerit.vecs_buff1, p->x, l, threading);
  apply_P_PRECISION(l->powerit.vecs_buff2, l->powerit.vecs_buff1, l, threading);
  vector_PRECISION_minus( l->powerit.vecs[i], p->x, l->powerit.vecs_buff2, start, end, l );
 }

 
 void blind_bp_op_PRECISION_apply( int op_id, level_struct* lx, struct Thread* threading ){

  //TODO: where to put this?
  // apply gamma5 before either operator
     int i;
  if( lx->powerit.spec_type == _SVs ){
    for( i=0;i<lx->powerit.nr_vecs;i++ ){
      if( lx->depth==0 ){
        gamma5_PRECISION( lx->powerit.vecs[i], lx->powerit.vecs[i], lx, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, lx->inner_vector_size, &startg5, &endg5, lx, threading, lx->num_lattice_site_var );
        coarse_gamma5_PRECISION( lx->powerit.vecs[i], lx->powerit.vecs[i], startg5, endg5, lx );
      }
    }
    SYNC_CORES(threading)
  }
  
  double buff_tol;

  complex_PRECISION* buff_b;
  complex_PRECISION* buff_x;
  gmres_PRECISION_struct* px = get_p_struct_PRECISION( lx );

  buff_tol = px->tol;
  buff_b = px->b;
  buff_x = px->x;

  START_MASTER(threading)
  px->tol = lx->powerit.bp_tol;
  END_MASTER(threading)

  double t0 = MPI_Wtime();
  for( int i=0;i<lx->powerit.nr_vecs;i++ ){

    START_MASTER(threading)
    px->b = lx->powerit.vecs[i];
    END_MASTER(threading)
    SYNC_CORES(threading)
      
    apply_to_one_vector(lx, i, threading);
    
    START_MASTER(threading)
    if (g.my_rank==0) printf(".");
    END_MASTER(threading)
  }
             // if(g.my_rank==0)printf("\n\nHiii--------------- %d of %d\n", 0, lx->powerit.nr_vecs);

  double t1 = MPI_Wtime();
  START_MASTER(threading)
  if (g.my_rank==0) printf("\n");
  //if (g.my_rank==0) printf("time block : %f\n", t1-t0);
  END_MASTER(threading)

  // restore values
  START_MASTER(threading)
  px->tol = buff_tol;
  px->b = buff_b;
  px->x = buff_x;
  END_MASTER(threading)

}



int apply_solver_powerit_PRECISION( level_struct* l, struct Thread *threading ){

    int nr_iters;
    double buff_coarsest_tol, buff_coarse_tol;
    
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    
    START_MASTER(threading)
    buff_coarse_tol = p->tol;
    p->tol = l->powerit.tol_buff;
    if( l->level==0 ){
      buff_coarsest_tol = g.coarse_tol;
      g.coarse_tol = l->powerit.tol_buff;
    }
    END_MASTER(threading)
    SYNC_CORES(threading)

    nr_iters = fgmres_PRECISION( p, l, threading );

    START_MASTER(threading)
    p->tol = buff_coarse_tol;
    if( l->level==0 ){
      g.coarse_tol = buff_coarsest_tol;
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)


    return nr_iters;
  }













 
        

        

        

