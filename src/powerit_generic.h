#ifndef POWERIT_PRECISION_HEADER
  #define POWERIT_PRECISION_HEADER

  void block_powerit_driver_PRECISION( level_struct* l, struct Thread* threading );
  void block_powerit_PRECISION_init_and_alloc( int spec_type, int op_id, int depth_bp_op, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading );
  void block_powerit_PRECISION_free( level_struct* l, struct Thread* threading );
  void block_powerit_PRECISION( int op_id, int depth_bp_op, level_struct *l, struct Thread* threading );
  void blind_bp_op_PRECISION_apply( op_id, lx, threading );
  void bp_qr_PRECISION( level_struct* lx, struct Thread* threading );
  void test_powerit_quality( int op_id, level_struct* lx, struct Thread* threading );


  void powerit_non_diff_op( level_struct *l, int i, struct Thread *threading );
  void powerit_diff_op( level_struct *l, int i, struct Thread *threading );
  void powerit_split_op( level_struct *l, int i, struct Thread *threading );

  int apply_solver_powerit_PRECISION( level_struct* l, struct Thread *threading );

  
  // this functional aplies anoperator to a bp vector 
  void (*apply_to_one_vector)();

#endif
