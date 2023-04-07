
    
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


  


  // -------------------------------------------------------------------

  
