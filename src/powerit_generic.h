#ifndef POWERIT_PRECISION_HEADER
  #define POWERIT_PRECISION_HEADER

  void block_powerit_driver_PRECISION( level_struct* l, struct Thread* threading );
  void block_powerit_PRECISION_init_and_alloc( int spec_type, int depth_bp_op_PRECISION, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading );
  void block_powerit_PRECISION_free( level_struct* l, struct Thread* threading );
  void block_powerit_PRECISION( int depth_bp_op_PRECISION, level_struct *l, struct Thread* threading );
  void blind_bp_op_PRECISION_apply( level_struct *lx, struct Thread *threading );
  void bp_qr_PRECISION( level_struct* lx, struct Thread* threading );
  void test_powerit_quality_PRECISION( level_struct* lx, struct Thread* threading );

  void powerit_non_diff_op_PRECISION( level_struct *l, int i, struct Thread *threading );
  void powerit_diff_op_PRECISION( level_struct *l, int i, struct Thread *threading );
  void powerit_split_full_rank_op_PRECISION( level_struct *l, int i, struct Thread *threading );
  void powerit_split_orthog_op_PRECISION( level_struct *l, int i, struct Thread *threading );

  int apply_solver_powerit_PRECISION( level_struct* l, struct Thread *threading );
  
  void get_rayleight_quotients_PRECISION(level_struct* l, struct Thread* threading );
  void compute_U_from_V_PRECISION( level_struct* l, struct Thread* threading );

#endif
