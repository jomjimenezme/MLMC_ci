#include "main.h"


// aux functions: from finest level, I get the level at desired depth
level_struct *get_correct_l_PRECISION( int depth_bp_op, level_struct* l ){
  int i;

  level_struct *lx = l;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = lx->next_level;
  }

  return lx;
}

// FIXME : hacked this function a bit to avoid compiler warnings
gmres_PRECISION_struct* get_p_struct_PRECISION_2( level_struct* l ){ 
  if( l->depth==0 ){
    if ( strcmp("PRECISION","float")==0 ) { return &(l->p_PRECISION); }
    else { return (gmres_PRECISION_struct*)&(g.p); }
  }
  else{ return &(l->p_PRECISION); }
}


// --------------------------------------------------------------------------------------------------------


void block_powerit_PRECISION_init_and_alloc( int spec_type, int depth_bp_op, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading ){
  int i;
  // access l at the right level
  level_struct* lx = get_correct_l_PRECISION( depth_bp_op,l );
  
  lx->powerit_PRECISION.nr_vecs = nr_vecs;
  lx->powerit_PRECISION.bp_tol = bp_tol;
  lx->powerit_PRECISION.nr_cycles = nr_bpi_cycles;

  lx->powerit_PRECISION.spec_type = spec_type;
  lx->powerit_PRECISION.gs_buffer = NULL;
  PUBLIC_MALLOC( lx->powerit_PRECISION.gs_buffer, complex_PRECISION, 2*lx->powerit_PRECISION.nr_vecs );

  lx->powerit_PRECISION.vecs = NULL;
  lx->powerit_PRECISION.U= NULL;
  PUBLIC_MALLOC( lx->powerit_PRECISION.vecs, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );
  PUBLIC_MALLOC( lx->powerit_PRECISION.U, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );
  lx->powerit_PRECISION.vecs[0] = NULL;
  lx->powerit_PRECISION.U[0] = NULL;
  PUBLIC_MALLOC( lx->powerit_PRECISION.vecs[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
  PUBLIC_MALLOC( lx->powerit_PRECISION.U[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    lx->powerit_PRECISION.vecs[i] = lx->powerit_PRECISION.vecs[0] + i*lx->vector_size;
    lx->powerit_PRECISION.U[i] = lx->powerit_PRECISION.U[0] + i*lx->vector_size;
  }

  lx->powerit_PRECISION.vecs_buff1 = NULL;
  PUBLIC_MALLOC( lx->powerit_PRECISION.vecs_buff1, complex_PRECISION, lx->vector_size );
  lx->powerit_PRECISION.vecs_buff2 = NULL;
  PUBLIC_MALLOC( lx->powerit_PRECISION.vecs_buff2, complex_PRECISION, lx->vector_size );
  lx->powerit_PRECISION.vecs_buff3 = NULL;
  PUBLIC_MALLOC( lx->powerit_PRECISION.vecs_buff3, complex_PRECISION, lx->vector_size );
  
  
  PUBLIC_MALLOC( lx->powerit_PRECISION.SV, complex_PRECISION, lx->powerit_PRECISION.nr_vecs  );
  PUBLIC_MALLOC( lx->powerit_PRECISION.EV, complex_PRECISION, lx->powerit_PRECISION.nr_vecs  ); 
  
  lx->powerit_PRECISION.M = malloc(lx->powerit_PRECISION.nr_vecs * lx->powerit_PRECISION.nr_vecs * sizeof(complex_PRECISION));
}


void block_powerit_PRECISION_free( level_struct* l, struct Thread* threading ){
  int j;

  for( j=0;j<g.num_levels;j++ ){
    // in case no deflation is requested
    if( g.trace_deflation_type[j]==3 ){ continue; }

    // access l at the right level
    level_struct* lx = get_correct_l_PRECISION( j,l );

    PUBLIC_FREE( lx->powerit_PRECISION.vecs[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
    PUBLIC_FREE( lx->powerit_PRECISION.vecs, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );
    PUBLIC_FREE( lx->powerit_PRECISION.U, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );

    PUBLIC_FREE( lx->powerit_PRECISION.vecs_buff1, complex_PRECISION, lx->vector_size );
    PUBLIC_FREE( lx->powerit_PRECISION.vecs_buff2, complex_PRECISION, lx->vector_size );
    PUBLIC_FREE( lx->powerit_PRECISION.vecs_buff3, complex_PRECISION, lx->vector_size );
    PUBLIC_FREE( lx->powerit_PRECISION.gs_buffer, complex_PRECISION, 2*lx->powerit_PRECISION.nr_vecs );
    
    PUBLIC_FREE( lx->powerit_PRECISION.SV, complex_PRECISION, lx->powerit_PRECISION.nr_vecs  );
    PUBLIC_FREE( lx->powerit_PRECISION.EV, complex_PRECISION, lx->powerit_PRECISION.nr_vecs  ); 
  
    free(lx->powerit_PRECISION.M);
  }
}


void block_powerit_driver_PRECISION( level_struct* l, struct Thread* threading ){
  int i, spec_type=0, depth_bp_op, nr_bp_vecs, nr_bpi_cycles;
  double bp_tol;
  level_struct *lx;

  // specify the following two params in the .ini input file, at different levels:
  // 	dx trace deflation type: 0   // 0 is difference, 1 is non-difference, 2 is split orthogonal, 3 is no deflation
  // 	dx trace deflation nr vectors

  printf("\n");
  for( i=0;i<g.num_levels;i++ ){
    if(g.my_rank==0)printf("BPI : level %d, request %d\n", i, g.trace_deflation_type[i]);

    depth_bp_op = i;
    nr_bp_vecs = g.trace_deflation_nr_vectors[i];
    bp_tol = g.trace_powerit_solver_tol[i];
    nr_bpi_cycles = g.trace_powerit_cycles[i];

    lx = get_correct_l_PRECISION( depth_bp_op,l );

    // in case no deflation is requested
    if( g.trace_deflation_type[i] == 0 ){ continue; }
    switch(g.trace_deflation_type[i]){
      case 1:
        lx->powerit_PRECISION.apply_to_one_vector = powerit_diff_op_PRECISION;
        break;
      case 2:
        lx->powerit_PRECISION.apply_to_one_vector = powerit_non_diff_op_PRECISION;
        break;
      case 3:
        lx->powerit_PRECISION.apply_to_one_vector = powerit_split_full_rank_op_PRECISION;
        break;
      case 4:
        lx->powerit_PRECISION.apply_to_one_vector = powerit_split_orthog_op_PRECISION;
        break;
      case 5:
        // big TODO !
        error0("under construction!\n");
        lx->powerit_PRECISION.apply_to_one_vector = powerit_split_orthog_op_PRECISION;
        break;
      case 6: //Deflation vetors for non-diff operator at the second level
        lx->next_level->powerit_PRECISION.apply_to_one_vector = powerit_non_diff_op_PRECISION;
        break;

      default:
        error0("Uknown type for operator in block power iteration\n");
    }

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
    //		   -- after calling power iteration, the result is in lx->powerit_PRECISION.vecs, with lx
    //		      the level struct of the chosen level
    if( depth_bp_op==(g.num_levels-1) && lx->powerit_PRECISION.apply_to_one_vector == powerit_diff_op_PRECISION){
      error0("There is no difference level operator at the coarsest level\n");
    }

    if(g.trace_deflation_type[i] != 6){ // same for all methods except multig. def.
        block_powerit_PRECISION_init_and_alloc( spec_type, depth_bp_op, nr_bp_vecs, nr_bpi_cycles, bp_tol, l, threading );
        block_powerit_PRECISION( depth_bp_op, l, threading );
    }else{ //Multigrid deflation
        //We need to allocate the memory on both finest and coarser level
        block_powerit_PRECISION_init_and_alloc( spec_type, depth_bp_op, nr_bp_vecs, nr_bpi_cycles, bp_tol, l, threading );
        block_powerit_PRECISION_init_and_alloc( spec_type, depth_bp_op+1, nr_bp_vecs, nr_bpi_cycles, bp_tol, l, threading );
        
        if(g.my_rank==0)
          printf("Doing BPI at depth %d \n", depth_bp_op+1);
        //and then call the method from the coarser level
        block_powerit_PRECISION( depth_bp_op+1, l, threading );
	    get_rayleight_quotients_PRECISION(depth_bp_op+1, l, threading);
        //test_orthogonality_PRECISION(depth_bp_op+1, l, threading );
        
	    //and compute the U vectors (also in coarser)
        compute_U_from_V_PRECISION( depth_bp_op+1, l, threading );
    }
    
  }
}


void matrix_computation_PRECISION(level_struct* lx, struct Thread* threading){
  int i, j;
  int k = lx->powerit_PRECISION.nr_vecs; 
  
  gmres_PRECISION_struct* p = get_p_struct_PRECISION_2( lx );
  
  //allocation of k buffer vectors
  vector_PRECISION* vec_buffer = NULL;
  PUBLIC_MALLOC( vec_buffer, complex_PRECISION*, k );
  vec_buffer[0] = NULL;
  PUBLIC_MALLOC( vec_buffer[0], complex_PRECISION, k*lx->vector_size );
  for( i=1; i<k; i++ ){
    vec_buffer[i] = vec_buffer[0] + i*lx->vector_size;
  }
  
  //vec_buffer = gamma_5 A V_c
  for( i=0; i<k; i++){
    apply_operator_PRECISION( vec_buffer[i], lx->powerit_PRECISION.vecs[i], p, lx, threading );
    if( lx->depth==0 ){
        gamma5_PRECISION(  vec_buffer[i],  vec_buffer[i], lx, threading );
    }
    else{
      int startg5, endg5;
      compute_core_start_end_custom(0, lx->inner_vector_size, &startg5, &endg5, lx, threading, lx->num_lattice_site_var );
      coarse_gamma5_PRECISION(  vec_buffer[i],  vec_buffer[i], startg5, endg5, lx );
    }
        
  }

  complex_PRECISION *M = lx->powerit_PRECISION.M;
  //M = U^* gamma_5 A V_c
  for( i=0; i<k; i++){
    for( j=0; j<k; j++){   
       M[i*k + j] = global_inner_product_PRECISION(lx->powerit_PRECISION.U[i], vec_buffer[j], p->v_start, p->v_end, lx, threading);
    }
  }
  /*//Example Matrix
    M[0] = 1;    M[1] = 0;    M[2] = 4;     // -5   0   -2
    M[3] = 1;    M[4] = 1;    M[5] = 6;     // -4   1   -1
    M[6] = -3;    M[7] = 0;    M[8] = -10;  // 3/2  0   1/2
  */
  // Call LAPACK function to invert the matrix M
  int ipiv[k]; // Pivot array
  //LU
  getrf_PRECISION(LAPACK_COL_MAJOR, k, k, M, k, ipiv);
  //Inverse
  getri_PRECISION(LAPACK_COL_MAJOR, k, M, k, ipiv);
  
  //Print matrix
   /* if(g.my_rank==0){
    printf("Inverted Matrix M:\n");
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
        printf("(%f) ", creal(M[i*k + j]));//, cimag(M[i*k + j]));
        }
        printf("\n");
    }
  }*/
  
  PUBLIC_FREE( vec_buffer, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );

}

void block_powerit_PRECISION( int depth_bp_op, level_struct *l, struct Thread *threading ){
  int i;

  // access l at the right level
  level_struct* lx = get_correct_l_PRECISION( depth_bp_op,l );

  // set the power iteration vectors to random
  START_LOCKED_MASTER(threading)
  vector_PRECISION_define_random( lx->powerit_PRECISION.vecs[0], 0, lx->powerit_PRECISION.nr_vecs*lx->vector_size, lx );
  END_LOCKED_MASTER(threading)
  SYNC_CORES(threading)

  for( i=0;i<lx->powerit_PRECISION.nr_cycles;i++ ){
    // apply the operator on the vectors ...
    blind_bp_op_PRECISION_apply( lx, threading );
    // ... and the resulting vectors are in lx->powerit_PRECISION.vecs
    //orthogonalize
    bp_qr_PRECISION( lx, threading );
  }

  // in the SVs case, this tests the eigenvectors coming out of the Hermitian problem
  //test_powerit_quality_PRECISION( op_id, lx, threading );
}


void bp_qr_PRECISION( level_struct* lx, struct Thread* threading ){
  // double call to Classical Gram Schmidt

  gram_schmidt_PRECISION( lx->powerit_PRECISION.vecs, lx->powerit_PRECISION.gs_buffer, 0, lx->powerit_PRECISION.nr_vecs, lx, threading );
  gram_schmidt_PRECISION( lx->powerit_PRECISION.vecs, lx->powerit_PRECISION.gs_buffer, 0, lx->powerit_PRECISION.nr_vecs, lx, threading );
}


// this is called in case we have chosen Plain Hutchinson
void powerit_non_diff_op_PRECISION( level_struct *l, int i, struct Thread *threading ){
  int start, end;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION_2( l );
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );

  apply_solver_powerit_PRECISION(l, threading );

  vector_PRECISION_copy( l->powerit_PRECISION.vecs[i], p->x, start, end, l );
}


// this is called in case we have chosen Traditional MGMLMC
void powerit_diff_op_PRECISION( level_struct *l, int i, struct Thread *threading ){
  int start, end;    

  gmres_PRECISION_struct* p = get_p_struct_PRECISION_2( l );
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );
    
  // coarse
  level_struct* lxc = l->next_level;
  gmres_PRECISION_struct* pxc = get_p_struct_PRECISION_2(lxc);
  apply_R_PRECISION(pxc->b, p->b, l, threading);     
  lxc->powerit_PRECISION.tol_buff = g.trace_powerit_solver_tol[lxc->depth];
  g.trace_powerit_solver_tol[lxc->depth] = g.trace_powerit_solver_tol[l->depth];
  apply_solver_powerit_PRECISION(lxc, threading);
  g.trace_powerit_solver_tol[lxc->depth] = lxc->powerit_PRECISION.tol_buff;
  apply_P_PRECISION(l->powerit_PRECISION.vecs_buff2, pxc->x, l, threading);
  
  // fine
  apply_solver_powerit_PRECISION(l, threading);
  vector_PRECISION_minus( l->powerit_PRECISION.vecs[i], p->x, l->powerit_PRECISION.vecs_buff2, start, end, l );
}


// big TODO : Generalize this
void powerit_split_orthog_op_PRECISION( level_struct *l, int i, struct Thread *threading ){
  int start, end;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );
  
 //Restrict. Prolongate, Difference, Solve, R, P, Difference.

  apply_R_PRECISION( l->powerit_PRECISION.vecs_buff1, l->powerit_PRECISION.vecs[i], l, threading );
  apply_P_PRECISION( l->powerit_PRECISION.vecs_buff2, l->powerit_PRECISION.vecs_buff1, l, threading );
  vector_PRECISION_minus(  p->b, l->powerit_PRECISION.vecs[i], l->powerit_PRECISION.vecs_buff2, start, end, l );
  apply_solver_powerit_PRECISION( l, threading );
  apply_R_PRECISION( l->powerit_PRECISION.vecs_buff1, p->x, l, threading );
  apply_P_PRECISION( l->powerit_PRECISION.vecs_buff2, l->powerit_PRECISION.vecs_buff1, l, threading );
  vector_PRECISION_minus(  l->powerit_PRECISION.vecs[i], p->x, l->powerit_PRECISION.vecs_buff2, start, end, l );
}


// the term tr( R A_{l}^{-1} P - A_{l+1}^{-1} )
// this is called in case we have chosen Split Full Rank MGMLMC
void powerit_split_full_rank_op_PRECISION( level_struct *l, int i, struct Thread *threading ){
  int start, end;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION_2( l );
  compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );

  apply_P_PRECISION( p->b, l->powerit_PRECISION.vecs[i], l, threading );
  apply_solver_powerit_PRECISION( l, threading );
  apply_R_PRECISION( l->powerit_PRECISION.vecs_buff2, p->x, l, threading );

  level_struct* lxc = l->next_level;
  gmres_PRECISION_struct* pxc = get_p_struct_PRECISION_2( lxc );
  compute_core_start_end( 0, lxc->inner_vector_size, &start, &end, lxc, threading );

  vector_PRECISION_copy( pxc->b, l->powerit_PRECISION.vecs[i], start, end, lxc );
  // solution of this solve is in l->next_level->p_PRECISION.x

  lxc->powerit_PRECISION.tol_buff = g.trace_powerit_solver_tol[lxc->depth];
  g.trace_powerit_solver_tol[lxc->depth] = g.trace_powerit_solver_tol[l->depth];
  apply_solver_powerit_PRECISION(lxc, threading);
  g.trace_powerit_solver_tol[lxc->depth] = lxc->powerit_PRECISION.tol_buff;

  vector_PRECISION_minus( l->powerit_PRECISION.vecs[i], l->powerit_PRECISION.vecs_buff2,  pxc->x,  start, end, l );
}


void blind_bp_op_PRECISION_apply( level_struct* lx, struct Thread* threading ){
  // apply gamma5 before any of the operators .. they're all Gamma5-Hermitian
  int i, start, end;

  // all the operators that we consider are Gamma5-Hermitian, so, we apply BPI on their
  // symmetrized versions. For this, then, we need to apply Gamma5 before the operator itself
  if( lx->powerit_PRECISION.spec_type == _SVs ){
    for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
      if( lx->depth==0 ){
        gamma5_PRECISION( lx->powerit_PRECISION.vecs[i], lx->powerit_PRECISION.vecs[i], lx, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, lx->inner_vector_size, &startg5, &endg5, lx, threading, lx->num_lattice_site_var );
        coarse_gamma5_PRECISION( lx->powerit_PRECISION.vecs[i], lx->powerit_PRECISION.vecs[i], startg5, endg5, lx );
      }
    }
    SYNC_CORES(threading)
  }

  gmres_PRECISION_struct* px = get_p_struct_PRECISION_2( lx );

  // after applying Gamma5, apply the operator on each deflation/BPI vector, sequentially
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    compute_core_start_end( 0, lx->inner_vector_size, &start, &end, lx, threading );
    vector_PRECISION_copy( px->b, lx->powerit_PRECISION.vecs[i], start, end, lx );

    lx->powerit_PRECISION.apply_to_one_vector(lx, i, threading);

    START_MASTER(threading)
    if (g.my_rank==0) { printf("."); fflush(0); }
    END_MASTER(threading)
  }

  START_MASTER(threading)
  if (g.my_rank==0) printf("\n");
  END_MASTER(threading)
}


int apply_solver_powerit_PRECISION( level_struct* l, struct Thread *threading ){
  int nr_iters;
  PRECISION buff_coarsest_tol=0, buff_coarse_tol=0;
    
  gmres_PRECISION_struct* p = get_p_struct_PRECISION_2( l );
    
  START_MASTER(threading)
  buff_coarse_tol = p->tol;
  p->tol = g.trace_powerit_solver_tol[l->depth];
  if( l->level==0 ){
    buff_coarsest_tol = g.coarse_tol;
    g.coarse_tol = g.trace_powerit_solver_tol[l->depth];
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

// U = \gamma_5 V Sign(Lambda)
void compute_U_from_V_PRECISION( int depth_bp_op, level_struct* l,  Thread* threading ){
  
  int i, start, end;
  level_struct* lx = get_correct_l_PRECISION( depth_bp_op,l );
  complex_PRECISION* lambda  = lx->powerit_PRECISION.EV;
  gmres_PRECISION_struct* px = get_p_struct_PRECISION_2( lx );
  compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

  PRECISION sign_lambda;
  for (i = 0; i< lx->powerit_PRECISION.nr_vecs; i++){
    if( lx->depth==0 ){
        gamma5_PRECISION( lx->powerit_PRECISION.U[i], lx->powerit_PRECISION.vecs[i], lx, threading );
    }
    else{
      int startg5, endg5;
      compute_core_start_end_custom(0, lx->inner_vector_size, &startg5, &endg5, lx, threading, lx->num_lattice_site_var );
      coarse_gamma5_PRECISION( lx->powerit_PRECISION.U[i], lx->powerit_PRECISION.vecs[i], startg5, endg5, lx );
    }
    
    if( signbit(creal(lambda[i])) ){
      sign_lambda = -1.0;
      vector_PRECISION_scale(lx->powerit_PRECISION.U[i], lx->powerit_PRECISION.U[i], sign_lambda, start, end, lx );
      //if(g.my_rank==0)
          //printf("Negative\n");
    }else{
     //if(g.my_rank==0)
          //printf("Positive\n");   
      continue;
    }

  }

}

void get_rayleight_quotients_PRECISION(int depth_bp_op, level_struct* l, struct Thread* threading ){
  
  int i, start, end;
  level_struct* lx = get_correct_l_PRECISION( depth_bp_op,l );
  complex_PRECISION rq = 0.0; PRECISION norm = 0.0;
  vector_PRECISION* vecs_buff;
  
  vecs_buff = NULL;
  PUBLIC_MALLOC( vecs_buff, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );
  vecs_buff[0] = NULL;  
  PUBLIC_MALLOC( vecs_buff[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
  
  START_MASTER(threading)
  for( i=1;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    vecs_buff[i] = vecs_buff[0] + i*lx->vector_size; //TODO:discuss inner_vector_size
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  
  gmres_PRECISION_struct* px = get_p_struct_PRECISION_2( lx );
  compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);
  
  // -----------------------------------
  // computing Rayleigh quotients

  // backup of deflation vectors
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    vector_PRECISION_copy( vecs_buff[i], lx->powerit_PRECISION.vecs[i], start, end, lx );
  }
  
  // apply the operator  (it includes Gamma5 for SV)
  blind_bp_op_PRECISION_apply( lx, threading );

  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    rq = global_inner_product_PRECISION( vecs_buff[i], lx->powerit_PRECISION.vecs[i],  px->v_start, px->v_end, lx, threading );
    norm = global_norm_PRECISION(vecs_buff[i], 0, lx->inner_vector_size, lx, threading );
    rq /= (norm*norm); 
    
    //We store the low modes of D:
    lx->powerit_PRECISION.EV[i] = 1.0/rq;
    lx->powerit_PRECISION.SV[i] = cabs_PRECISION(1.0/rq);
    // Restore Backup
    vector_PRECISION_copy(lx->powerit_PRECISION.vecs[i], vecs_buff[i], start, end, lx );
    
    if(g.my_rank==0)
        printf("---------\t %f + i%f\n", creal(rq), cimag(rq));
  }
  
  // -----------------------------------
  // computing eigen-residuals
  
  PRECISION norm2, rel_res;

  // backup of deflation vectors
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    vector_PRECISION_copy( vecs_buff[i], lx->powerit_PRECISION.vecs[i], start, end, lx );
  }
  
  // apply the operator
  blind_bp_op_PRECISION_apply( lx, threading );
 
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    rq = 1.0/lx->powerit_PRECISION.EV[i];
    
    // Ax - \lambda x
    vector_PRECISION_saxpy( lx->powerit_PRECISION.vecs[i], lx->powerit_PRECISION.vecs[i],
                            vecs_buff[i], -rq, start, end, lx );// z := x + alpha*y

      
    norm  = global_norm_PRECISION(lx->powerit_PRECISION.vecs[i], 0, lx->inner_vector_size, lx, threading );
    norm2 = global_norm_PRECISION(vecs_buff[i], 0, lx->inner_vector_size, lx, threading );
    
    rel_res = norm/(norm2*cabs_PRECISION(rq)); // ||Ax - \lambda x||/|| \lambda x ||
    
    
    
    // Restore Backup
    vector_PRECISION_copy(lx->powerit_PRECISION.vecs[i], vecs_buff[i], start, end, lx );
    
    if(g.my_rank==0)
        printf("---------\t %f\n", rel_res);
  }
  
  PUBLIC_FREE( vecs_buff[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
}




void test_orthogonality_PRECISION(int depth_bp_op, level_struct* l, struct Thread* threading ){
  
  int i, j, start, end;
  level_struct* lx = get_correct_l_PRECISION( depth_bp_op,l );
  complex_PRECISION aux;

  gmres_PRECISION_struct* px = get_p_struct_PRECISION_2( lx );
  compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);
  
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    for( j =0; j< lx->powerit_PRECISION.nr_vecs; j++){
      aux = global_inner_product_PRECISION( lx->powerit_PRECISION.vecs[i], lx->powerit_PRECISION.vecs[j], start, end, lx, threading    );
      if(g.my_rank==0) printf(" %.12f  ", creal(aux));
    }
      if(g.my_rank==0) printf("\n");
  }
  
}
// ------------------------------------------------------------------------------------------


void test_powerit_quality_PRECISION( level_struct* lx, struct Thread* threading ){
  int i, start, end;
  vector_PRECISION* vecs_buff1;
  vector_PRECISION* vecs_buff2;

  gmres_PRECISION_struct* px = get_p_struct_PRECISION_2( lx );

  compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

  vecs_buff1 = NULL;
  vecs_buff2 = NULL;
  PUBLIC_MALLOC( vecs_buff2, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );
  PUBLIC_MALLOC( vecs_buff1, complex_PRECISION*, lx->powerit_PRECISION.nr_vecs );
  
  vecs_buff1[0] = NULL;
  vecs_buff2[0] = NULL;  
  PUBLIC_MALLOC( vecs_buff1[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
  PUBLIC_MALLOC( vecs_buff2[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
  
  START_MASTER(threading)
  for( i=1;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    vecs_buff1[i] = vecs_buff1[0] + i*lx->vector_size;
    vecs_buff2[i] = vecs_buff2[0] + i*lx->vector_size;
  }
  END_MASTER(threading)
  SYNC_CORES(threading)

  // backup of lx->powerit_PRECISION.vecs
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    vector_PRECISION_copy( vecs_buff1[i], lx->powerit_PRECISION.vecs[i], start, end, lx );
  }

  // apply the operator
  blind_bp_op_PRECISION_apply( lx, threading );

  // swap pointers  
  complex_PRECISION** buff_ptr = lx->powerit_PRECISION.vecs;
  complex_PRECISION** buff_ptr1 = vecs_buff1;    
  START_MASTER(threading)
  lx->powerit_PRECISION.vecs = buff_ptr1;
  vecs_buff1 = buff_ptr;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    // compute the Rayleigh quotient
    complex_PRECISION rq;
    rq = global_inner_product_PRECISION( lx->powerit_PRECISION.vecs[i], vecs_buff1[i], px->v_start, px->v_end, lx, threading );
    double norm = global_norm_PRECISION( lx->powerit_PRECISION.vecs[i], 0, lx->inner_vector_size, lx, threading );
    rq /= norm;

    // print the Rayleigh quotient
    START_MASTER(threading)
    if (g.my_rank==0) printf( "Rayleigh quotient = %.16f+i%.16f \t",CSPLIT(rq) );
    END_MASTER(threading)

    // compute the eigenvalue residual
    vector_PRECISION_scale( vecs_buff2[i], lx->powerit_PRECISION.vecs[i], rq, start, end, lx );
    vector_PRECISION_minus( vecs_buff1[i], vecs_buff1[i], vecs_buff2[i], start, end, lx );
    double resx = global_norm_PRECISION( vecs_buff1[i], 0, lx->inner_vector_size, lx, threading ) / ( cabs_PRECISION(rq)*norm );
    
    // print the residuals
    START_MASTER(threading)
    if (g.my_rank==0) printf( "Eigenvalue residual = %.16f\n",resx );
    END_MASTER(threading)
  }
  
  //Restore backup
  for( i=0;i<lx->powerit_PRECISION.nr_vecs;i++ ){
    vector_PRECISION_copy( lx->powerit_PRECISION.vecs[i], vecs_buff1[i], start, end, lx );
  }
  //swap pointers Back
  START_MASTER(threading)
  vecs_buff1 = buff_ptr1;
  lx->powerit_PRECISION.vecs = buff_ptr;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  PUBLIC_FREE( vecs_buff1[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
  PUBLIC_FREE( vecs_buff2[0], complex_PRECISION, lx->powerit_PRECISION.nr_vecs*lx->vector_size );
}
