
    
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
        
        //MLMC
        h->mlmc_b1 = NULL;
        
        h->mlmc_b2 =NULL;

        h->mlmc_testing =NULL;

        
        h->rademacher_vector =NULL;
        
        SYNC_MASTER_TO_ALL(threading)
    }
    void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading ) {
        hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;
        
        //For MLMC
        PUBLIC_MALLOC( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );
        PUBLIC_MALLOC( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );
        PUBLIC_MALLOC( h->mlmc_testing, complex_PRECISION, l->inner_vector_size );
        
        
        PUBLIC_MALLOC( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );
        
    }
    
    void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading ) {
        hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;
        
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
    complex_PRECISION one_sample, variance, trace;
    double RMSD;
    struct sample estimate;

    // TODO : move this allocation to some init function
    complex_PRECISION* samples = (complex_PRECISION*) malloc( h->max_iters*sizeof(complex_PRECISION) );
    memset( samples, 0.0, h->max_iters*sizeof(complex_PRECISION) );

    estimate.acc_trace = 0.0;

    double t0 = MPI_Wtime();
    for( i=0; i<h->max_iters;i++ ){

      // 1. create Rademacher vector, stored in h->rademacher_vector
      rademacher_create_PRECISION( l, h, type, threading );

      // 2. apply the operator to the Rademacher vector
      // 3. dot product
      one_sample = h->hutch_compute_one_sample( l, h, threading );
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
        if(g.my_rank==0) printf( "%d\tVariance:\t%f\n", i, creal(variance));
        END_MASTER(threading);
        RMSD = sqrt(creal(variance)/j);
        if( i > h->min_iters && RMSD < cabs(trace) * h->trace_tol * h->tol_per_level[l->depth]) break; 
      }
    }
    double t1 = MPI_Wtime();

    START_MASTER(threading);


    if(g.my_rank==0) printf( "%d\t \tvariance = %f+i%f \t t = %f, \t d = %.3f\n", i, CSPLIT(variance), t1-t0, h->tol_per_level[l->depth]);


    END_MASTER(threading);

    estimate.sample_size = i;

    free(samples);

    return estimate;
  }


  gmres_PRECISION_struct* get_p_struct_PRECISION( level_struct* l ){

    if( l->depth==0 ){
      return &(g.p);
    }
    else{
      return &(l->p_PRECISION);
    }
  }

  int apply_solver_PRECISION( level_struct* l, struct Thread *threading ){

    int nr_iters;
    double buff1, buff2;

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
    
    double t0 = MPI_Wtime();
    nr_iters = fgmres_PRECISION( p, l, threading );
    double t1 = MPI_Wtime();

   /* int start, end;
    compute_core_start_end( p->v_start, p->v_end, &start, &end, l, threading );
    apply_operator_PRECISION( p->w, p->x, p, l, threading ); // compute w = D*x
    vector_PRECISION_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
    double norm = global_norm_PRECISION( p->r, p->v_start, p->v_end, l, threading ); //||r||
    double relative = norm/global_norm_PRECISION( p->b, p->v_start, p->v_end, l, threading );//||r||/b

    if(g.my_rank==0)printf("-----------------------------------\n-----------------------------------\n");
	if(g.my_rank==0)printf("\t Solve time %f,\t Iters %d, \t ||r||= %e, \t ||r||/||b|| = %e, \t used_tol %e,\t coarsest:tol %e\n", t1-t0, nr_iters, norm, relative, p->tol, g.coarse_tol);
    if(g.my_rank==0)printf("-----------------------------------\n-----------------------------------\n");
*/

    if(g.my_rank==0)printf("\t Solve time %f,\t Iters %d\n", t1-t0, nr_iters);
    START_MASTER(threading);
    p->tol = buff1;
    if( l->level==0 ){
      g.coarse_tol = buff2;
    }
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);

    return nr_iters;
  }


  complex_PRECISION hutchinson_plain_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  
    {
      int start, end;
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    }
    double t0 = MPI_Wtime();
    {
      apply_solver_PRECISION( l, threading );
    }
    double t1 = MPI_Wtime();
    START_MASTER(threading);
      if(g.my_rank==0) printf( "Plain Solve Time:\t%f\n", t1-t0);
    END_MASTER(threading);
    // subtract the results and perform dot product
    {
      int start, end;
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      complex_double aux = global_inner_product_double( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );   
          if(g.my_rank==0)  printf( "\t----> one-level solve <-----\t%f \n", creal(aux) );
        return aux;  
    }
  }


  complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "Trace computation via Hutchinson's method ...\n" );
    END_MASTER(thrading);
    
    int i;
    complex_PRECISION trace = 0.0;
    struct sample estimate;
    hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
    level_struct* lx;
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tfinest (and only) level ...\n" );
    END_MASTER(thrading);

    lx = l;
    // set the pointer to the finest-level Hutchinson estimator
    h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
    
    double t_plain = MPI_Wtime();

    estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
    
    trace += estimate.acc_trace/estimate.sample_size;
    double t_plain1 = MPI_Wtime();
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    if(g.my_rank==0)  printf( "Average solve time %f \n", (t_plain1-t_plain)/h->max_iters );
    END_MASTER(thrading);

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "... done\n" );
    END_MASTER(thrading);

    return trace;
  }



  // apply the interpolation
  void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){

    if( l->depth==0 ){
      interpolate3_PRECISION( l->sbuf_PRECISION[0], in, l, threading );
      trans_back_PRECISION( out, l->sbuf_PRECISION[0], l->s_PRECISION.op.translation_table, l, threading );
    }
    else{
      interpolate3_PRECISION( out, in, l, threading );
    }
  }


  // apply the restriction
  void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){

    if( l->depth==0 ){
      trans_PRECISION( l->sbuf_PRECISION[0], in, l->s_PRECISION.op.translation_table, l, threading );     
      restrict_PRECISION( out, l->sbuf_PRECISION[0], l, threading );
    }
    else{
      restrict_PRECISION( out, in, l, threading );
    }
  }


 // the term tr( A_{l}^{-1} - P A_{l+1}^{-1} R )
  complex_PRECISION hutchinson_mlmc_difference_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

    // FIRST TERM : result stored in p->x
    // apply A_{l}^{-1}
    {
      int start, end;
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
      // solution of this solve is in l->p_PRECISION.x
      apply_solver_PRECISION( l, threading );
    }

    // SECOND TERM : result stored in h->mlmc_b2
    // 1. Restrict
    // 2. invert
    // 3. Prolongate
    {
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
      
      apply_R_PRECISION( l->next_level->p_PRECISION.b, h->rademacher_vector, l, threading );
      // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
      apply_solver_PRECISION( l->next_level, threading );
      apply_P_PRECISION( h->mlmc_b2, l->next_level->p_PRECISION.x, l, threading );
    }

    // subtract the results and perform dot product
    {
      int start, end;
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l);
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_PRECISION_minus( h->mlmc_b1, p->x, h->mlmc_b2, start, end, l); 

      //Deflate here?
      complex_double aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l, threading );
          if(g.my_rank==0)  printf( "\t----> Difference-level solve <-----\t%f \n", creal(aux) );
        return aux; 
    }
  }



  complex_PRECISION mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "Trace computation via 'traditional' difference levels ...\n" );
    END_MASTER(thrading);

    int i;
    complex_PRECISION trace = 0.0;
    struct sample estimate;
    hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
    level_struct* lx;

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tdifference levels ...\n" );
    END_MASTER(thrading);
    
    // for all but coarsest level
    lx = l->next_level->next_level;
    for( i=2; i<g.num_levels-1;i++ ){
      // set the pointer to the mlmc difference operator
      h->hutch_compute_one_sample = hutchinson_mlmc_difference_PRECISION;
      estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
      trace += estimate.acc_trace/estimate.sample_size;
      lx = lx->next_level;
    }
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tcoarsest level ...\n" );
    END_MASTER(thrading);
    
    // coarsest level
    // set the pointer to the coarsest-level Hutchinson estimator
    h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
    estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "... done\n" );
    END_MASTER(thrading);

    return trace;
  }


 complex_PRECISION split_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "Trace computation via split levels ...\n" );
    END_MASTER(thrading);

    int i;
    complex_PRECISION trace = 0.0;
    struct sample estimate;
    hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
    level_struct* lx;
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tfull-rank difference levels ...\n" );
    END_MASTER(thrading);
    
    // for all but coarsest level
    lx = l;
    for( i=0; i<g.num_levels-1;i++ ){  
      // set the pointer to the split intermediate operator
      h->hutch_compute_one_sample = hutchinson_split_intermediate_PRECISION;
      estimate = hutchinson_blind_PRECISION( lx, h, 1, threading );
      trace += estimate.acc_trace/estimate.sample_size;
      lx = lx->next_level;    
    }
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);

    START_MASTER(threading);
    if(g.my_rank==0) printf( "\torthogonalized difference levels ...\n" );
    END_MASTER(thrading);
    
    // for all but coarsest level
    lx = l;
    for( i=0; i<g.num_levels-1;i++ ){      
      // set the pointer to the split orthogonal operator
      h->hutch_compute_one_sample = hutchinson_split_orthogonal_PRECISION;
      estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
      trace += estimate.acc_trace/estimate.sample_size;
      lx = lx->next_level; 
    }
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tcoarsest level ...\n" );
    END_MASTER(thrading);

    // coarsest level
    // set the pointer to the coarsest-level Hutchinson estimator
    h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
    estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);

    return trace;
  }


  // the term tr( A_{l}^{-1}(I - P_{l} P_{l}^{H})  )
  complex_PRECISION hutchinson_split_orthogonal_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

    // 1. project
    // 2. invert

    // FIRST TERM

    {
      int start, end;
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

      apply_R_PRECISION( h->mlmc_b2, h->rademacher_vector, l, threading );
      apply_P_PRECISION( h->mlmc_b1, h->mlmc_b2, l, threading );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_PRECISION_minus( h->mlmc_b1, h->rademacher_vector, h->mlmc_b1, start, end, l );

      vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l );
    }
    
    // SECOND "factor"

    {
      apply_solver_PRECISION( l, threading );
    }

    // perform dot product
    {
      int start, end;
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      
      complex_double aux = global_inner_product_PRECISION( h->mlmc_b1, p->x, p->v_start, p->v_end, l, threading );        
      if(g.my_rank==0)  printf( "\t----> Orthogonal-level solve <-----\t%f \n", creal(aux) );
      return aux; 
      
    }
  }



 // the term tr( R A_{l}^{-1} P - A_{l+1}^{-1} )
  complex_PRECISION hutchinson_split_intermediate_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

    // FIRST TERM : result stored in h->mlmc_b1

    // 1. prolongate
    // 2. invert
    // 3. restrict
    {
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

      apply_P_PRECISION( p->b, h->rademacher_vector, l, threading );
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
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l->next_level );
      // solution of this solve is in l->next_level->p_PRECISION.x
      apply_solver_PRECISION( l->next_level, threading );
      vector_PRECISION_copy( h->mlmc_b2, l->next_level->p_PRECISION.x, start, end, l->next_level );
    }

    // subtract the results and perform dot product
    {
      int start, end;
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l->next_level );
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start, &end, l->next_level, threading );
      vector_PRECISION_minus( h->mlmc_b1, h->mlmc_b1, h->mlmc_b2, start, end, l->next_level ); 
      
      complex_double aux = global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l->next_level, threading );         
      if(g.my_rank==0)  printf( "\t----> Intermediate-level solve <-----\t%f \n", creal(aux) );
      return aux; 
      //return global_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l->next_level, threading );   
    }
    
  }
  


  // -------------------------------------------------------------------

  
