#include "main.h"

struct estimate {
  int counter; //required number of estimates.
  complex_PRECISION estimate;
};

struct sample {
  // required number of estimates
  int sample_size;
  // accumulated trace
  complex_PRECISION acc_trace;
};

    
void hutchinson_diver_PRECISION_init( level_struct *l, struct Thread *threading ) {
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);

  h->max_iters = NULL;
  h->min_iters = NULL;
        
  // MLMC
  h->mlmc_b1 = NULL;
  h->mlmc_b2 =NULL;
  h->mlmc_testing =NULL;
  h->rademacher_vector =NULL;
  SYNC_MASTER_TO_ALL(threading)
}


void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading ) {
  int i;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;

  PUBLIC_MALLOC( h->max_iters, int, g.num_levels );
  PUBLIC_MALLOC( h->min_iters, int, g.num_levels );

  //For MLMC
  PUBLIC_MALLOC( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->mlmc_testing, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );

  for ( i=0;i<g.num_levels;i++ ) {
    h->max_iters[i] = g.trace_max_iters[i];
    h->min_iters[i] = g.trace_min_iters[i];
  }
}


void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading ) {
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;

  PUBLIC_FREE( h->max_iters, int, g.num_levels );
  PUBLIC_FREE( h->min_iters, int, g.num_levels );
        
  PUBLIC_FREE( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );   
  PUBLIC_FREE( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );   
  PUBLIC_FREE( h->mlmc_testing, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );   
}


void rademacher_create_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading ){
  if( type==0 ){
    START_MASTER(threading)
    vector_PRECISION_define_random_rademacher( h->rademacher_vector, 0, l->inner_vector_size, l );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
  }
  else if( type==1 ){
    START_MASTER(threading)
    vector_PRECISION_define_random_rademacher( h->rademacher_vector, 0, l->next_level->inner_vector_size, l->next_level );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
   }
  else{ error0("Unknown value for type of Rademacher vector in relation to level of creation\n"); }
}


// type : in case of 0 create Rademacher vectors at level l, in case of 1 create Rademacher vectors at level l->next_level
struct sample hutchinson_blind_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading ){
  int i, j;
  complex_PRECISION one_sample=0, variance=0, trace=0;
  double RMSD;
  struct sample estimate;

  // TODO : move this allocation to some init function
  complex_PRECISION* samples = (complex_PRECISION*) malloc( h->max_iters[l->depth]*sizeof(complex_PRECISION) );
  memset( samples, 0.0, h->max_iters[l->depth]*sizeof(complex_PRECISION) );

  estimate.acc_trace = 0.0;

  for( i=0; i<h->max_iters[l->depth];i++ ){
    // 1. create Rademacher vector, stored in h->rademacher_vector
    rademacher_create_PRECISION( l, h, type, threading );

    // 2. apply the operator to the Rademacher vector
    // 3. dot product
    one_sample = h->hutch_compute_one_sample( -1, l, h, threading );

    samples[i] = one_sample;

    // 4. compute estimated trace and variance, print something?
    estimate.acc_trace += one_sample;

    if( i!=0 ){
      variance = 0.0;
      estimate.sample_size = i+1;
      trace = estimate.acc_trace/estimate.sample_size;
      for( j=0; j<i; j++ ){
        variance += conj(samples[j] - trace) * (samples[j] - trace);
      }
      variance = variance / j;
      START_MASTER(threading);
      if(g.my_rank==0) {
        printf("[%d, trace: %f+%f, variance: %f] ", i, creal(trace), cimag(trace), creal(variance));
        fflush(0);
      }
      END_MASTER(threading);
      RMSD = sqrt(creal(variance)/j);
      if( i > h->min_iters[l->depth] && RMSD < cabs(trace) * h->trace_tol * h->tol_per_level[l->depth]) break; 
    }
  }
  if(g.my_rank==0) printf("\n");

  estimate.sample_size = i;

  free(samples);

  return estimate;
}


// FIXME : hacked this function a bit to avoid compiler warnings
gmres_PRECISION_struct* get_p_struct_PRECISION( level_struct* l ){
  if( l->depth==0 ){
    if ( strcmp("PRECISION","float")==0 ) { return &(l->p_PRECISION); }
    else { return (gmres_PRECISION_struct*)&(g.p); }
  }
  else{ return &(l->p_PRECISION); }
}


int apply_solver_PRECISION( level_struct* l, struct Thread *threading ){
  int nr_iters;
  double buff1=0, buff2=0;

  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    
  buff1 = p->tol;
  p->tol = g.tol;
  START_MASTER(threading);
  if( l->level==0 ){
    buff2 = g.coarse_tol;
    g.coarse_tol = g.tol;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);
    
  nr_iters = fgmres_PRECISION( p, l, threading );

  START_MASTER(threading);
  p->tol = buff1;
  if( l->level==0 ){
    g.coarse_tol = buff2;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);

  return nr_iters;
}


// if type_appl==-1 : Hutchinson-like
// else : direct term, where type_appl is the index of the deflation vector to apply
//        the operator on
complex_PRECISION hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    } else {
      vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }
  }

  {
    apply_solver_PRECISION( l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      if(g.trace_deflation_type[l->depth] != 0){
        hutchinson_deflate_vector_PRECISION(p->x, l, threading); 
      }
      aux = global_inner_product_PRECISION( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );   
    } else {
      vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);  
      aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    }

    return aux;  
  }
}


complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){    
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
  trace += estimate.acc_trace/estimate.sample_size;

  // if deflation vectors are available
  if(g.trace_deflation_type[l->depth] != 0){
    trace += hutchinson_deflated_direct_term_PRECISION( l, h, threading );
  }

  return trace;
}


// apply the interpolation
void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){
  if( l->depth==0 ){
    interpolate3_PRECISION( l->sbuf_PRECISION[0], in, l, threading );
    trans_back_PRECISION( (vector_double)out, l->sbuf_PRECISION[0], l->s_PRECISION.op.translation_table, l, threading );
  }
  else{
    interpolate3_PRECISION( out, in, l, threading );
  }
}


// apply the restriction
void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){
  if( l->depth==0 ){
    trans_PRECISION( l->sbuf_PRECISION[0], (vector_double)in, l->s_PRECISION.op.translation_table, l, threading );     
    restrict_PRECISION( out, l->sbuf_PRECISION[0], l, threading );
  }
  else{
    restrict_PRECISION( out, in, l, threading );
  }
}


// the term tr( A_{l}^{-1} - P A_{l+1}^{-1} R )
complex_PRECISION hutchinson_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // FIRST TERM : result stored in p->x
  // apply A_{l}^{-1}
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    } else {
      vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // solution of this solve is in l->p_PRECISION.x
    apply_solver_PRECISION( l, threading );
  }

  // SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict
  // 2. invert
  // 3. Prolongate
  {

    if ( type_appl==-1 ) {
      apply_R_PRECISION( l->next_level->p_PRECISION.b, h->rademacher_vector, l, threading );
    } else {
      apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );
    apply_P_PRECISION( h->mlmc_b2, l->next_level->p_PRECISION.x, l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l);
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_minus( h->mlmc_b1, p->x, h->mlmc_b2, start, end, l); 

    if ( type_appl==-1 ) {
      if(g.trace_deflation_type[l->depth] != 0){
        hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading); 
      }
      aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l, threading );
    } else {
      aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }

    return aux; 
  }
}


complex_PRECISION mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1; i++ ){ 
    // set the pointer to the mlmc difference operator
    h->hutch_compute_one_sample = hutchinson_mlmc_difference_PRECISION;
    estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;
    // if deflation vectors are available
    if(g.trace_deflation_type[lx->depth] != 0){
      trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
    }
    lx = lx->next_level;
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
  trace += estimate.acc_trace/estimate.sample_size;
  // if deflation vectors are available
  if(g.trace_deflation_type[lx->depth] != 0){
    trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
  }

  return trace;
}


complex_PRECISION split_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx=0;

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1 ;i++ ){  
    // set the pointer to the split full rank operator
    h->hutch_compute_one_sample = hutchinson_split_intermediate_PRECISION;
    estimate = hutchinson_blind_PRECISION( lx, h, 1, threading );
    trace += estimate.acc_trace/estimate.sample_size;

    // if deflation vectors are available
    if( g.trace_deflation_type[lx->depth] != 0 ){
      if( g.trace_deflation_type[lx->depth]==3 || g.trace_deflation_type[lx->depth]==5 ){
        trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
      }
    }
    lx = lx->next_level;    
  }

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1;i++ ){      
    // set the pointer to the split orthogonal operator
    h->hutch_compute_one_sample = hutchinson_split_orthogonal_PRECISION;
    estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;

    // if deflation vectors are available
    if( g.trace_deflation_type[lx->depth] != 0 ){
      if( g.trace_deflation_type[lx->depth]==4 || g.trace_deflation_type[lx->depth]==5 ){
        trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
      }
    }
    lx = lx->next_level; 
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
  // if deflation vectors are available
  if(g.trace_deflation_type[lx->depth] != 0){
    trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
  }
  trace += estimate.acc_trace/estimate.sample_size;

  return trace;
}


// the term tr( A_{l}^{-1}(I - P_{l} P_{l}^{H})  )
complex_PRECISION hutchinson_split_orthogonal_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // 1. project
  // 2. invert

  // FIRST TERM
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

    if ( type_appl==-1 ) {
      apply_R_PRECISION( h->mlmc_b2, h->rademacher_vector, l, threading );
    } else {
      apply_R_PRECISION( h->mlmc_b2, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    apply_P_PRECISION( h->mlmc_b1, h->mlmc_b2, l, threading );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_minus( p->b, h->rademacher_vector, h->mlmc_b1, start, end, l );
    } else {
      vector_PRECISION_minus( p->b, l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, start, end, l );
    }
  }
    
  // SECOND "factor"
  {
    apply_solver_PRECISION( l, threading );
  }

  // perform dot product
  {
    int start, end;
    complex_PRECISION aux = 0.0;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      
    
    apply_R_PRECISION( h->mlmc_b2, p->x, l, threading );
    apply_P_PRECISION( h->mlmc_b1, h->mlmc_b2, l, threading );

    vector_PRECISION_minus( h->mlmc_b1, p->x, h->mlmc_b1, start, end, l );
    
    if ( type_appl==-1 ) {
      if(g.trace_deflation_type[l->depth] != 0){
        hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading); 
      }
      aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l, threading );   
    } else {  
      aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
    return aux; 
  }
}


// the term tr( R A_{l}^{-1} P - A_{l+1}^{-1} )
complex_PRECISION hutchinson_split_intermediate_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

  // FIRST TERM : result stored in h->mlmc_b1
  // 1. prolongate
  // 2. invert
  // 3. restrict
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

    if ( type_appl==-1 ) {
      apply_P_PRECISION( p->b, h->rademacher_vector, l, threading );
    } else {
      apply_P_PRECISION( p->b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }
    // the input of this solve is p->x, the output p->b
    apply_solver_PRECISION( l, threading );
    apply_R_PRECISION( h->mlmc_b1, p->x, l, threading );
  }

  // SECOND TERM : result stored in h->mlmc_b2

  // apply A_{l+1}^{-1}
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l->next_level );
    compute_core_start_end( 0, l->next_level->inner_vector_size, &start, &end, l->next_level, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l->next_level );
    } else {
      vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l->next_level );
    }

    // solution of this solve is in l->next_level->p_PRECISION.x
    apply_solver_PRECISION( l->next_level, threading );
    vector_PRECISION_copy( h->mlmc_b2, l->next_level->p_PRECISION.x, start, end, l->next_level );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l->next_level );
    compute_core_start_end( 0, l->next_level->inner_vector_size, &start, &end, l->next_level, threading );
    vector_PRECISION_minus( h->mlmc_b1, h->mlmc_b1, h->mlmc_b2, start, end, l->next_level ); 

    if(g.trace_deflation_type[l->depth] != 0){
      hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading); 
    }

    if ( type_appl==-1 ) {
      aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l->next_level, threading );         
    } else {
      aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l->next_level, threading );         
    }
      
    return aux;
  }
}


complex_PRECISION hutchinson_deflated_direct_term_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  int i, start, end;
  complex_PRECISION direct_trace = 0.0;

  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
  // 0. create small matrix to store all dot products, let's call it small_T
  complex_PRECISION small_T[l->powerit_PRECISION.nr_vecs];

  for( i=0; i<l->powerit_PRECISION.nr_vecs;i++ ){
    // 1. apply the operator on the ith deflation vector
    small_T[i] = h->hutch_compute_one_sample( i, l, h, threading );

    // 3. take trace of small_T, store in estimate.direct_trace
    direct_trace += small_T[i];
  }
  if(g.my_rank==0) printf("Trace from the direct term : %f +i %f\n", CSPLIT(direct_trace));
  return direct_trace;
}


void hutchinson_deflate_vector_PRECISION(vector_PRECISION input, level_struct *l, struct Thread *threading ){
  int i, start, end;
  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

  complex_PRECISION aux[l->powerit_PRECISION.nr_vecs];

  for( i=0; i<l->powerit_PRECISION.nr_vecs; i++ ){
    aux[i] = global_inner_product_PRECISION(l->powerit_PRECISION.vecs[i], input, p->v_start, p->v_end, l, threading);	
  }

  // buff1 <- alpha0 * v0
  vector_PRECISION_scale( l->powerit_PRECISION.vecs_buff1 , l->powerit_PRECISION.vecs[0], aux[0], start, end, l);
  for( int i = 1;  i< l->powerit_PRECISION.nr_vecs; i++ ){
    // buff3 <- buff1
    vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff3, l->powerit_PRECISION.vecs_buff1, start, end, l);
    // buff2 <- alphai * vi
    vector_PRECISION_scale( l->powerit_PRECISION.vecs_buff2, l->powerit_PRECISION.vecs[i], aux[i], start, end, l);
    // buff1 <- buff3 + buff2
    vector_PRECISION_plus( l->powerit_PRECISION.vecs_buff1 , l->powerit_PRECISION.vecs_buff3 , l->powerit_PRECISION.vecs_buff2, start, end, l);
  }

  vector_PRECISION_minus(  input, input, l->powerit_PRECISION.vecs_buff1, start, end, l );
}




complex_PRECISION multigrid_deflation_driver_PRECISION( level_struct *l, struct Thread *threading ){


  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;
  level_struct* lxc;
  
  lx = l;
  lxc = lx->next_level;
  //Precompute k-by-k matrix 
  matrix_computation_PRECISION( lxc,  threading);
  
  // set the pointer to the deflated operator
  h->hutch_compute_one_sample = hutchinson_multigrid_deflated_PRECISION;
  estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
  trace += estimate.acc_trace/estimate.sample_size;

  //adding direct term
  if(g.trace_deflation_type[lx->depth] == 6){
    trace += hutchinson_multigrid_direct_PRECISION(l, h, threading );
  }
  
  if(g.my_rank==0) printf("TRACE: %f+%f\n", CSPLIT(trace));
  return trace;
}

  


complex_PRECISION hutchinson_multigrid_deflated_PRECISION(int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
    
   
    int i, j, start, end, startxc, endxc, k;
    complex_PRECISION aux =0.0 +0.0*I;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION_2( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    level_struct* lxc = l->next_level;
    compute_core_start_end( 0, lxc->inner_vector_size, &startxc, &endxc, lxc, threading );
   
    k = l->powerit_PRECISION.nr_vecs;
    complex_PRECISION *M_inv = lxc->powerit_PRECISION.M;
    complex_PRECISION* vec_buffer = NULL;
    vec_buffer = malloc(lxc->powerit_PRECISION.nr_vecs * sizeof(complex_PRECISION));
    complex_PRECISION *vec_buffer1 = NULL;
    vec_buffer1 = malloc(k * sizeof(complex_PRECISION));
    memset(vec_buffer, 0.0 + 0.0*I, k * sizeof(complex_PRECISION));
    memset(vec_buffer1, 0.0 + 0.0*I, k * sizeof(complex_PRECISION));
    
    //gamma_5 x
    if( l->depth==0 ){
        gamma5_PRECISION(  l->powerit_PRECISION.vecs_buff1,  h->rademacher_vector, l, threading );
    }
    else{
      int startg5, endg5;
      compute_core_start_end_custom(0, l->inner_vector_size, &startg5, &endg5, l, threading, l->num_lattice_site_var );
      coarse_gamma5_PRECISION(  l->powerit_PRECISION.vecs_buff1,  h->rademacher_vector, startg5, endg5, l );
    }
    
    //R gamma_5 x
    apply_R_PRECISION(lxc->powerit_PRECISION.vecs_buff1, l->powerit_PRECISION.vecs_buff1, l, threading);
     //U_c R gamma_5 x
    for( i=0; i<k; i++){   
       vec_buffer[i] = global_inner_product_PRECISION(lxc->powerit_PRECISION.U[i], lxc->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, lxc, threading);
    }
    
   
   //vec_buffer1 = M_inv vec_buffer = M_inv U_c R gamma_5 x
   for (i = 0; i < k; ++i) {
        for (j = 0; j < k; ++j) {
            vec_buffer1[i] += M_inv[i * k + j] * vec_buffer[j]; // Compute y[i] += A[i][j] * x[j]
        }
    }
  
   //V_c M_inv U_c R gamma_5 x
    vector_PRECISION_define( lxc->powerit_PRECISION.vecs_buff1, 0.0 + 0.0*I, startxc, endxc, lxc );
    for( i = 0; i < lxc->powerit_PRECISION.nr_vecs; i++ ){
    vector_PRECISION_saxpy( lxc->powerit_PRECISION.vecs_buff1 , lxc->powerit_PRECISION.vecs_buff1,
                            lxc->powerit_PRECISION.vecs[i], vec_buffer1[i], startxc, endxc, lxc );
    }
   
   //P V_c M_inv U_c R gamma_5 x
   apply_P_PRECISION(l->powerit_PRECISION.vecs_buff1, lxc->powerit_PRECISION.vecs_buff1, l, threading);
    
        
   vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    

    // solution of this solve is in l->p_PRECISION.x
    apply_solver_PRECISION( l, threading );
  
    vector_PRECISION_minus( l->powerit_PRECISION.vecs_buff1, p->x, l->powerit_PRECISION.vecs_buff1, start, end, l );
    
    aux = global_inner_product_PRECISION( h->rademacher_vector, l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    
    return aux;

 
    
}




//M_inv (U_c \gamma_5c V_c)
complex_PRECISION hutchinson_multigrid_direct_PRECISION(level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

  int i, j, startx, endx, startg5, endg5,  k;
  //everything is done in coarser level
  level_struct* lxc = l->next_level;
  gmres_PRECISION_struct* pxc = get_p_struct_PRECISION_2( lxc );
  compute_core_start_end(0, lxc->inner_vector_size, &startx, &endx, lxc, threading);
  compute_core_start_end_custom(0, lxc->inner_vector_size, &startg5, &endg5, lxc, threading, lxc->num_lattice_site_var );
  k = lxc->powerit_PRECISION.nr_vecs;
  
  //Buffer
  vector_PRECISION* vecs_buff;
  vecs_buff = NULL;
  PUBLIC_MALLOC( vecs_buff, complex_PRECISION*, lxc->powerit_PRECISION.nr_vecs );
  vecs_buff[0] = NULL;
  PUBLIC_MALLOC( vecs_buff[0], complex_PRECISION, lxc->powerit_PRECISION.nr_vecs*lxc->vector_size );

  START_MASTER(threading)
  for( i=1; i<lxc->powerit_PRECISION.nr_vecs; i++ ){
    vecs_buff[i] = vecs_buff[0] + i*lxc->vector_size; //TODO:discuss inner_vector_size
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
 
 
  // copy of deflation vectors
  for( i=0; i<lxc->powerit_PRECISION.nr_vecs; i++ ){
    vector_PRECISION_copy( vecs_buff[i], lxc->powerit_PRECISION.vecs[i], startx, endx, lxc );
  }
  
  //gamma_5c V_c
  for(i=0; i<k; i++){
    coarse_gamma5_PRECISION(  vecs_buff[i],  lxc->powerit_PRECISION.vecs[i], startg5, endg5, lxc );
  }
 
  complex_PRECISION *C = NULL;
  C = malloc(k * k * sizeof(complex_PRECISION));
 
  for(i=0; i<k; i++){
    for(j=0; j<k; j++){
      C[i + j * k] = global_inner_product_PRECISION( lxc->powerit_PRECISION.U[i], vecs_buff[j],  pxc->v_start, pxc->v_end, lxc, threading );
    }
  }

  complex_PRECISION *M_inv = lxc->powerit_PRECISION.M; 
  complex_PRECISION *aux = malloc(k * k * sizeof(complex_PRECISION));
   

  matrix_multiply_PRECISION(aux, M_inv, C, k);  
  
  complex_PRECISION trace = 0.0 + 0.0*I; 
  
  for(i=0; i<k*k; i+=k+1){
    trace += aux[i];
  }
  
  if(g.my_rank==0) printf("Trace from direct term = %f +%f \n", CSPLIT(trace));
  return trace; 

}

//helper function: temporary till working with BLAS C=A*B
void matrix_multiply_PRECISION(complex_PRECISION *C, complex_PRECISION *A, complex_PRECISION *B, int n) {
/*if(g.my_rank==0)printf("Matrix A:\n");
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          if(g.my_rank==0)printf("(%.2f, %.2f)\t", creal(A[ i*n +j ]), cimag(A[i*n + j]));
      }
      if(g.my_rank==0)printf("\n");
  }

if(g.my_rank==0)printf("B Matrix:\n");
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          if(g.my_rank==0)printf("(%.2f, %.2f)\t", creal(B[ i*n +j ]), cimag(B[i*n + j]));
      }
      if(g.my_rank==0)printf("\n");
  }
*/

for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        C[i * n + j] = 0.0 + 0.0*I;
        for (int k = 0; k < n; k++) {
            C[i * n + j] += A[i * n + k] * B[j + k*n]; // Note the change here
        }
    }
} 
     

/*
if(g.my_rank==0)printf("Matrix A(B) :\n");
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          if(g.my_rank==0)printf("(%.2f, %.2f)\t", creal(C[ i*n +j ]), cimag(C[i*n + j]));
      }
      if(g.my_rank==0)printf("\n");
  }

*/

}
