#ifndef HUTCHINSON_PRECISION_HEADER
  #define HUTCHINSON_PRECISION_HEADER

  struct Thread;

  complex_PRECISION mlmc_hutchinson_diver_PRECISION( level_struct *l, struct Thread *threading );

  void hutchinson_diver_PRECISION_init( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading );

  void rademacher_create_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading );
   
  complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION split_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  
  void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );
  void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );
  int apply_solver_powerit_PRECISION( level_struct* l, struct Thread *threading );

  complex_PRECISION hutchinson_deflated_direct_term_PRECISION(level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading);
  void hutchinson_deflate_vector_PRECISION(vector_PRECISION input, level_struct *l, struct Thread *threading );

  complex_PRECISION hutchinson_deflated_direct_term_difference_PRECISION(level_struct *l, struct Thread *threading);

  // different operators for Hutchinson and block power iteration
  complex_PRECISION hutchinson_split_orthogonal_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  complex_PRECISION hutchinson_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  complex_PRECISION hutchinson_split_intermediate_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );

  complex_PRECISION hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
 
  int apply_solver_PRECISION( level_struct* l, struct Thread *threading );
  gmres_PRECISION_struct* get_p_struct_PRECISION( level_struct* l );
  struct sample hutchinson_blind_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type, struct Thread *threading );
  
  complex_PRECISION hutchinson_multigrid_deflated_PRECISION(int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  complex_PRECISION hutchinson_multigrid_direct_PRECISION(level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  void matrix_multiply_PRECISION(complex_PRECISION *C, complex_PRECISION *A, complex_PRECISION *B, int n);
#endif
