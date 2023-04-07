
    
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
        
        h->rademacher_vector =NULL;
        
        SYNC_MASTER_TO_ALL(threading)
    }
    void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading ) {
        hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;
        
        //For MLMC
        PUBLIC_MALLOC( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );
        PUBLIC_MALLOC( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );
        
        PUBLIC_MALLOC( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );
        
    }
    
    void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading ) {
        hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;
        
        PUBLIC_FREE( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );   
        PUBLIC_FREE( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );   
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
      return global_inner_product_PRECISION( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );   
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




  // -------------------------------------------------------------------

  
