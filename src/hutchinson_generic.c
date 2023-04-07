
    
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
    

    
    



  // -------------------------------------------------------------------

  
