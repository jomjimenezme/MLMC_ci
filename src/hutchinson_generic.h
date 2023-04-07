#ifndef HUTCHINSON_PRECISION_HEADER
  #define HUTCHINSON_PRECISION_HEADER


  struct Thread;

  complex_PRECISION mlmc_hutchinson_diver_PRECISION( level_struct *l, struct Thread *threading );

  void hutchinson_diver_PRECISION_init( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading );
  
  complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION split_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );

  // different operators for Hutchinson and block power iteration
  complex_PRECISION hutchinson_split_orthogonal_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  
  
  void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );
  void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );

  int apply_solver_powerit_PRECISION( level_struct* l, struct Thread *threading );
#endif
