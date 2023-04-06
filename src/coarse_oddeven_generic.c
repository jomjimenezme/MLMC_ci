/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#include "main.h"

void coarse_selfcoupling_LU_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l ) {

  // clover = [ A B      , A=A*, D=D*, C = -B*
  //            C D ]
  //
  // order: upper triangle of A, upper triangle of D, B, each column major
  //
  // tm_term = [ E 0      , E=-E*, F=-F* diag. excluded
  //             0 F ]
  //
  // order: upper triangle of E, upper triangle of F
  //
  // output = [ A+E  B   
  //             C  D+F ] LU decomposed

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n;
  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n2+1)*j] = *clover;
    clover++; // diagonal entry
  }
  // D
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n2+1)*j] = *clover;
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = -conj_PRECISION(*clover);
      clover++;      
    }
  }

#ifdef HAVE_TM
  config_PRECISION tm_term = op->tm_term + n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n2*i+j] += *tm_term;
        output[i+n2*j] += -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n2+1)*j] += *tm_term;
      tm_term++; // diagonal entry
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
        output[n2*i+j] += *tm_term;
        output[i+n2*j] += -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n2+1)*j] += *tm_term;
      tm_term++; // diagonal entry
    }
  }
#endif
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n2; k++ ) {
    for ( i=k+1; i<n2; i++ ) {
      output[n2*i+k] = output[n2*i+k]/output[(n2+1)*k]; // output(i,k) = output(i,k)/output(k,k)
      for ( j=k+1; j<n2; j++ )
        output[n2*i+j] = output[n2*i+j]-output[n2*i+k]*output[n2*k+j]; // output(i,j) = output(i,j)-output(i,k)*output(k,j)
    }
  }
}

#ifdef HAVE_TM1p1
void coarse_selfcoupling_LU_doublet_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l ) {

  // clover = [ A B      , A=A*, D=D*, C = -B*
  //            C D ]
  //
  // order: upper triangle of A, upper triangle of D, B, each column major
  //
  // tm_term = [ E 0      , E=-E*, F=-F* diag. excluded
  //             0 F ]
  //
  // order: upper triangle of E, upper triangle of F
  //
  // epsbar_term = [ G 0      , G=-G*, H=-H* diag. excluded
  //                 0 H ]
  //
  // order: upper triangle of G, upper triangle of H
  //
  // output = [ A+E  G   B   0
  //             G  A-E  0   B
  //             C   0  D+F  H
  //             0   C   H  D-F ]  LU decomposed

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n, n3 = 3*n, n4 = 4*n;
  // set the matrix up
  // 0
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*(i+0 )+(j+n3)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n )+(j+n2)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n2)+(j+n )] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n3)+(j+0 )] = _COMPLEX_PRECISION_ZERO;
    }
  }

  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // D
  for ( j=n2; j<n3; j++ ) {
    for ( i=n2; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n2; j<n3; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = -conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = -conj_PRECISION(*clover);
      clover++;      
    }
  }

#ifdef HAVE_TM
  config_PRECISION tm_term = op->tm_term + n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n4*i+j] += *tm_term;
        output[i+n4*j] += -conj_PRECISION(*tm_term);
        output[n4*(i+n)+(j+n)] -= *tm_term;
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n4+1)*j] += *tm_term;
      output[(n4+1)*(j+n)] -= *tm_term;
      tm_term++; // diagonal entry
    }
    // F
    for ( j=n2; j<n3; j++ ) {
      for ( i=n2; i<j; i++ ) {
        output[n4*i+j] += *tm_term;
        output[i+n4*j] += -conj_PRECISION(*tm_term);
        output[n4*(i+n)+(j+n)] -= *tm_term;
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n4+1)*j] += *tm_term;
      output[(n4+1)*(j+n)] -= *tm_term;
      tm_term++; // diagonal entry
    }
  }
#endif

  config_PRECISION epsbar_term = op->epsbar_term + n*(n+1)*index;
  // G
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
  // H
  for ( j=n3; j<n4; j++ ) {
    for ( i=n2; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n4; k++ ) {
    for ( i=k+1; i<n4; i++ ) {
      output[n4*i+k] = output[n4*i+k]/output[(n4+1)*k]; // output(i,k) = output(i,k)/output(k,k)
      for ( j=k+1; j<n4; j++ )
        output[n4*i+j] = output[n4*i+j]-output[n4*i+k]*output[n4*k+j]; // output(i,j) = output(i,j)-output(i,k)*output(k,j)
    }
  }
}
#endif


void coarse_perform_fwd_bwd_subs_PRECISION( vector_PRECISION x, vector_PRECISION b, config_PRECISION A, level_struct *l ) {
  
  register int i, j, n2 = l->num_lattice_site_var;
  
  // solve x = U^(-1) L^(-1) b
  // forward substitution with L
  for ( i=0; i<n2; i++ ) {
    x[i] = b[i];
    for ( j=0; j<i; j++ ) {
      x[i] = x[i] - A[i*n2+j]*x[j];
    }
  }
  // backward substitution with U
  for ( i=n2-1; i>=0; i-- ) {
    for ( j=i+1; j<n2; j++ ) {
      x[i] = x[i] - A[i*n2+j]*x[j];
    }
    x[i] = x[i]/A[i*(n2+1)];
  }
}


void coarse_LU_multiply_PRECISION( vector_PRECISION y, vector_PRECISION x, config_PRECISION A, level_struct *l ) {
  
  register int i, j, n2 = l->num_lattice_site_var;
  
  // y = Ax
  // multiplication with U
  for ( i=0; i<n2; i++ ) {
    y[i] = A[i*(n2+1)]*x[i];
    for ( j=i+1; j<n2; j++ )
      y[i] += A[i*n2+j]*x[j];
  }
  // multiplication with L
  for ( i=n2-1; i>0; i-- )
    for ( j=0; j<i; j++ )
      y[i] += A[i*n2+j]*y[j];
}


void coarse_diag_ee_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start, end;
  compute_core_start_end_custom( 0, op->num_even_sites, &start, &end, l, threading, 1 );
  // even sites
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION
  coarse_self_couplings_PRECISION( y, x, op, start, end, l );
#else
  coarse_self_couplings_PRECISION_vectorized( y, x, op, start, end, l );
#endif
}

void coarse_diag_oo_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start, end;
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION 
  int num_site_var=l->num_lattice_site_var,
    oo_inv_size = SQUARE(num_site_var);
#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else
  config_PRECISION sc = op->clover_oo_inv;
#endif

  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );

  x += num_site_var*(op->num_even_sites+start);
  y += num_site_var*(op->num_even_sites+start);  
  sc += oo_inv_size*start;

  for ( int i=start; i<end; i++ ) {
    coarse_LU_multiply_PRECISION( y, x, sc, l );
    x += num_site_var;
    y += num_site_var;
    sc += oo_inv_size;
  }
  
#else
  compute_core_start_end_custom( op->num_even_sites, l->num_inner_lattice_sites, &start, &end, l, threading, 1 );
  coarse_self_couplings_PRECISION_vectorized( y, x, op, start, end, l );
#endif
}

void coarse_diag_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l ) {
  
  coarse_diag_ee_PRECISION( y, x, op, l, no_threading );
  coarse_diag_oo_PRECISION( y, x, op, l, no_threading );
}

void coarse_diag_oo_inv_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, 
                               level_struct *l, struct Thread *threading ) {
  
  int start, end;
  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );
  
  // odd sites
  int num_site_var = l->num_lattice_site_var,
    oo_inv_size = SQUARE(num_site_var);

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION
#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else
  config_PRECISION sc = op->clover_oo_inv;
#endif
#else
  int lda = SIMD_LENGTH_PRECISION*((num_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  oo_inv_size = 2*num_site_var*lda;
#ifdef HAVE_TM1p1
  OPERATOR_TYPE_PRECISION *sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv_vectorized:op->clover_oo_inv_vectorized;
#else
  OPERATOR_TYPE_PRECISION *sc = op->clover_oo_inv_vectorized;
#endif
#endif

  x += num_site_var*(op->num_even_sites+start);
  y += num_site_var*(op->num_even_sites+start);  
  sc += oo_inv_size*start;

  for ( int i=start; i<end; i++ ) {
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION
    coarse_perform_fwd_bwd_subs_PRECISION( y, x, sc, l );
#else
    for(int j=0; j<num_site_var; j++)
      y[j] = _COMPLEX_PRECISION_ZERO;
    cgemv( num_site_var, sc, lda, (float *)x, (float *)y);
#endif
    x += num_site_var;
    y += num_site_var;
    sc += oo_inv_size;
  }
}


void coarse_oddeven_PRECISION_set_self_couplings( level_struct *l, struct Thread *threading ) {

  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
  int nv = l->num_parent_eig_vect, start, end;

  coarse_operator_PRECISION_set_self_couplings( op, l, threading );
  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1);

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION

  int size = SQUARE(2*nv);
  for( int i=start; i<end; i++ )
    coarse_selfcoupling_LU_decomposition_PRECISION( op->clover_oo_inv+i*size, op, op->num_even_sites+i, l );

#ifdef HAVE_TM1p1
  int size_doublet = SQUARE(4*nv);
  for( int i=start; i<end; i++ )
    coarse_selfcoupling_LU_doublet_decomposition_PRECISION( op->clover_doublet_oo_inv+i*size_doublet, op, 
                                                            op->num_even_sites+i, l );
#endif

// this else is for vectorized
#else

  int column_offset = SIMD_LENGTH_PRECISION*((2*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int size_v = 2*2*nv*column_offset;
  for( int i=start; i<end; i++ )
    cgem_inverse( 2*nv, op->clover_oo_inv_vectorized + i*size_v, 
                  op->clover_vectorized + (op->num_even_sites+i)*size_v, column_offset );

#ifdef BLOCK_JACOBI
  // preparation for direct solves in Block Jacobi
  {
    // l->p_PRECISION.block_jacobi_PRECISION.bj_op_inv_vectorized
    // l->p_PRECISION.block_jacobi_PRECISION.bj_op_vectorized

    // 1. prepare the operator : B = D_{ee} - I_{eo}*D_{oo}^{-1}*I_{oe}
    OPERATOR_TYPE_PRECISION* Dee = op->clover_vectorized;
    OPERATOR_TYPE_PRECISION* Dooinv = op->clover_oo_inv_vectorized;
    int startx,endx;
    compute_core_start_end_custom( 0, op->num_even_sites, &startx, &endx, l, threading, 1);
    // size of bj_op_vectorized : 2*2*nv*column_offset*op->num_even_sites
    PRECISION sc_fctr = 0.0;
    for( int i=startx; i<endx; i++ ){
      for( int ix=0;ix<size_v;ix++ )
        l->p_PRECISION.block_jacobi_PRECISION.bj_op_vectorized[i*size_v+ix] = Dee[i*size_v+ix] - sc_fctr*Dooinv[i*size_v+ix];
    }

    // 2. invert B
    for( int i=startx; i<endx; i++ ){
      cgem_inverse( 2*nv, l->p_PRECISION.block_jacobi_PRECISION.bj_op_inv_vectorized + i*size_v, 
                    l->p_PRECISION.block_jacobi_PRECISION.bj_op_vectorized + i*size_v, column_offset );
    }
  }
#endif

#ifdef HAVE_TM1p1
  int column_doublet_offset = SIMD_LENGTH_PRECISION*((4*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int size_doublet_v = 2*4*nv*column_doublet_offset;
  for( int i=start; i<end; i++ )
    cgem_inverse( 4*nv, op->clover_doublet_oo_inv_vectorized + i*size_doublet_v, 
                  op->clover_doublet_vectorized + (op->num_even_sites+i)*size_doublet_v, column_doublet_offset );

#ifdef BLOCK_JACOBI
  // preparation for direct solves in Block Jacobi
  {
    // l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_inv_vectorized
    // l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_vectorized

    // 1. prepare the operator : B = D_{ee} - I_{eo}*D_{oo}^{-1}*I_{oe}
    OPERATOR_TYPE_PRECISION* Dee = op->clover_doublet_vectorized;
    OPERATOR_TYPE_PRECISION* Dooinv = op->clover_doublet_oo_inv_vectorized;
    int startx,endx;
    compute_core_start_end_custom( 0, op->num_even_sites, &startx, &endx, l, threading, 1);
    // size of bj_op_vectorized : 2*2*nv*column_offset*op->num_even_sites
    PRECISION sc_fctr = 0.0;
    for( int i=startx; i<endx; i++ ){
      for( int ix=0;ix<size_doublet_v;ix++ )
        l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_vectorized[i*size_doublet_v+ix] = Dee[i*size_doublet_v+ix] - sc_fctr*Dooinv[i*size_doublet_v+ix];
    }

    // 2. invert B
    for( int i=startx; i<endx; i++ ){
      cgem_inverse( 4*nv, l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_inv_vectorized + i*size_doublet_v, 
                    l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_vectorized + i*size_doublet_v, column_doublet_offset );
    }
  }
#endif

#endif
#endif
}

void coarse_oddeven_PRECISION_set_couplings( level_struct *l, struct Thread *threading ) {

  coarse_oddeven_PRECISION_set_self_couplings( l, threading );
  coarse_operator_PRECISION_set_neighbor_couplings( &(l->oe_op_PRECISION), l, threading );

}

void coarse_oddeven_alloc_PRECISION( level_struct *l ) {

  int nv = l->num_parent_eig_vect, oe_offset=0, mu, **bt = NULL,
    *eot = NULL, *nt = NULL, *tt = NULL, t, z, y, x, le[4], N[4];
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  operator_PRECISION_alloc( op, _ODDEVEN, l );

  // buffers
  MALLOC( op->buffer, complex_PRECISION*, 2 );
  op->buffer[0] = NULL;
#ifdef HAVE_TM1p1
  MALLOC( op->buffer[0], complex_PRECISION, 4*l->vector_size );
  op->buffer[1] = op->buffer[0] + 2*l->vector_size;  
#else
  MALLOC( op->buffer[0], complex_PRECISION, 2*l->vector_size );
  op->buffer[1] = op->buffer[0] + l->vector_size;  
#endif

  for ( mu=0; mu<4; mu++ ) {
    le[mu] = l->local_lattice[mu];
    N[mu] = le[mu]+1;
    op->table_dim[mu] = N[mu];
  }

  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;

  // estimate site numbers
  op->num_even_sites = 0;
  op->num_odd_sites = 0;
  op->oe_offset = oe_offset;
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            op->num_odd_sites++;
          } else {
            op->num_even_sites++;
          }
        }
  
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION

  MALLOC( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#ifdef HAVE_TM1p1
  MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
#endif

#else
  int column_offset = SIMD_LENGTH_PRECISION*((2*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  MALLOC_HUGEPAGES( op->clover_oo_inv_vectorized, PRECISION, 2*2*nv*column_offset*op->num_odd_sites, 4*SIMD_LENGTH_PRECISION );
#ifdef BLOCK_JACOBI
  l->p_PRECISION.block_jacobi_PRECISION.bj_op_inv_vectorized = NULL;
  MALLOC_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_op_inv_vectorized, PRECISION, 2*2*nv*column_offset*op->num_even_sites, 4*SIMD_LENGTH_PRECISION );
  l->p_PRECISION.block_jacobi_PRECISION.bj_op_vectorized = NULL;
  MALLOC_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_op_vectorized, PRECISION, 2*2*nv*column_offset*op->num_even_sites, 4*SIMD_LENGTH_PRECISION );
#endif
#ifdef HAVE_TM1p1
  int column_doublet_offset = SIMD_LENGTH_PRECISION*((4*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  MALLOC_HUGEPAGES( op->clover_doublet_oo_inv_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_odd_sites, 4*SIMD_LENGTH_PRECISION );
#ifdef BLOCK_JACOBI
  l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_inv_vectorized = NULL;
  MALLOC_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_inv_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_even_sites, 4*SIMD_LENGTH_PRECISION );
  l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_vectorized = NULL;
  MALLOC_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_even_sites, 4*SIMD_LENGTH_PRECISION );
#endif
#endif

#endif

  // define data layout
  eot = op->index_table;
  define_eot( eot, N, l );

  // neighbor table, translation table
  nt = op->neighbor_table;
  tt = op->translation_table;
  define_nt_bt_tt( nt, op->backward_neighbor_table, NULL, tt, eot, N, l );

  // boundary table
  bt = op->c.boundary_table;
  define_eo_bt( bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites, op->c.num_boundary_sites, N, l );

  // ghost
  ghost_sendrecv_init_PRECISION( _COARSE_GLOBAL, &(op->c), l ) ;

  // solver
  if ( l->level == 0 )
    l->p_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
  else
    l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;

#ifdef BLOCK_JACOBI
  // solver
  if ( l->level == 0 )
    l->p_PRECISION.block_jacobi_PRECISION.local_p.v_end = op->num_even_sites*l->num_lattice_site_var;
  //else
  //  l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
#endif

}

void coarse_oddeven_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l, 
                                     struct Thread *threading ) {

  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  START_LOCKED_MASTER(threading)
    int ns=l->num_inner_lattice_sites, nv = l->num_parent_eig_vect, i,
    D_size = 4*SQUARE(2*nv),
    clover_size = (nv)*(nv*2+1),
    block_size = (nv)*(nv+1);
  config_PRECISION D_in = in->D,
    clover_in = in->clover,
    odd_proj_in = in->odd_proj;

  // neighbor couplings
  if ( reorder ) {
    int t, z, y, x, index, *le = l->local_lattice, oe_offset = op->oe_offset,
      *it = in->index_table, *dt = in->table_dim;
    config_PRECISION D_oe = op->D, 
      D_eo = (op->D)+D_size*op->num_even_sites,
      clover_ee = op->clover,
      clover_oo = (op->clover)+clover_size*op->num_even_sites,
      odd_proj_ee = op->odd_proj,
      odd_proj_oo = op->odd_proj+block_size*op->num_even_sites;

    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            index = site_index( t, z, y, x, dt, it );
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              for ( i=0; i<D_size; i++ ) 
                D_eo[i] = D_in[ index*D_size+i ];
              for ( i=0; i<clover_size; i++ )
                clover_oo[i] = clover_in[ index*clover_size+i ];
              for ( i=0; i<block_size; i++ )
                odd_proj_oo[i] = odd_proj_in[ index*block_size+i ];
              D_eo += D_size;
              clover_oo += clover_size;
              odd_proj_oo += block_size;
            } else {
              for ( i=0; i<D_size; i++ )
                D_oe[i] = D_in[ index*D_size+i ];
              for ( i=0; i<clover_size; i++ )
                clover_ee[i] = clover_in[ index*clover_size+i ];
              for ( i=0; i<block_size; i++ )
                odd_proj_ee[i] = odd_proj_in[ index*block_size+i ];
              D_oe += D_size;
              clover_ee += clover_size;
              odd_proj_ee += block_size;
            }
          }
    
  } else {
    for ( i=0; i<D_size*ns; i++ )
      op->D[i] = D_in[i];
    for ( i=0; i<clover_size*ns; i++ )
      op->clover[i] = clover_in[i];
    for ( i=0; i<block_size*ns; i++ ) {
      op->odd_proj[i] = odd_proj_in[i];
    }
    
  }
  END_LOCKED_MASTER(threading)
  
  op->m0 = in->m0;

#ifdef HAVE_TM
  tm_term_PRECISION_setup( in->mu, in->mu_even_shift, in->mu_odd_shift, op, l, threading );
#endif  
#ifdef HAVE_TM1p1
  epsbar_term_PRECISION_setup( in->epsbar, in->epsbar_ig5_even_shift, in->epsbar_ig5_odd_shift, op, l, threading );
#endif
  
  coarse_oddeven_PRECISION_set_couplings( l, threading );
  
}


void coarse_oddeven_free_PRECISION( level_struct *l ) {
  
  int nv = l->num_parent_eig_vect, vs = l->vector_size;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  operator_PRECISION_free( op, _ODDEVEN, l );
  coarse_operator_PRECISION_free_vectorized( op, l );

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION

  FREE( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#ifdef HAVE_TM1p1
  FREE( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
#endif

#else
  int column_offset = SIMD_LENGTH_PRECISION*((2*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  FREE_HUGEPAGES( op->clover_oo_inv_vectorized, PRECISION, 2*2*nv*column_offset*op->num_odd_sites );
#ifdef BLOCK_JACOBI
  FREE_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_op_inv_vectorized, PRECISION, 2*2*nv*column_offset*op->num_even_sites );
  FREE_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_op_vectorized, PRECISION, 2*2*nv*column_offset*op->num_even_sites );
#endif
#ifdef HAVE_TM1p1
  int column_doublet_offset = SIMD_LENGTH_PRECISION*((4*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  FREE_HUGEPAGES( op->clover_doublet_oo_inv_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_odd_sites );
#ifdef BLOCK_JACOBI
  FREE_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_inv_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_even_sites );
  FREE_HUGEPAGES( l->p_PRECISION.block_jacobi_PRECISION.bj_doublet_op_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_even_sites );
#endif
#endif

#endif
  
#ifdef HAVE_TM1p1
  FREE( op->buffer[0], complex_PRECISION, 4*vs );
#else
  FREE( op->buffer[0], complex_PRECISION, 2*vs );
#endif
  FREE( op->buffer, complex_PRECISION*, 2 );
}


void coarse_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                    const int amount, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      num_4link_var=4*4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      num_link_var=4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

  int core_start;
  int core_end;

#ifdef PERS_COMMS
  g.pers_comms_id1 = 0;
#endif
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );
  
  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  //START_MASTER(threading)
  //printf0("print 2 ...\n");
  //END_MASTER(threading)

  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  //START_MASTER(threading)
  //printf0("print 3 ...\n");
  //END_MASTER(threading)

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}


void coarse_n_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading ) {

#ifdef PERS_COMMS
  g.pers_comms_id1 = 1;
#endif

#ifdef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
#ifndef COMM_HIDING_COARSEOP
  int sign = -1;
  coarse_pn_hopping_term_PRECISION_vectorized( out, in, op, amount, l, sign, threading);
#else
  coarse_n_hopping_term_PRECISION_vectorized( out, in, op, amount, l, threading );
#endif
  return;
#else
  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      num_4link_var=4*4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      num_link_var=4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );
  
  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  
  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                    const int amount, level_struct *l, struct Thread *threading ) {

#ifdef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;

  OPERATOR_TYPE_PRECISION *D_vectorized;
  int column_offset = 2*SIMD_LENGTH_PRECISION*((l->num_parent_eig_vect+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int vectorized_link_offset = 2*2*l->num_parent_eig_vect*column_offset;

  int core_start;
  int core_end;

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );

  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)

  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 0*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 1*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 2*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 3*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_vectorized + 4*vectorized_link_offset*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    D_vectorized += vectorized_link_offset;
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    in_pt = in + num_site_var*op->neighbor_table[index+X];
    D_vectorized += vectorized_link_offset;
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_pn_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                    const int amount, level_struct *l, int sign, struct Thread *threading ) {

#ifdef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION

  START_NO_HYPERTHREADS(threading)

  int mu, i, num_site_var=l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;

  OPERATOR_TYPE_PRECISION *D_vectorized;
  int column_offset = 2*SIMD_LENGTH_PRECISION*((l->num_parent_eig_vect+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int link_offset = 2*2*l->num_parent_eig_vect*column_offset;
  int *neighbor_fw = op->neighbor_table;
  int *neighbor_bw = op->backward_neighbor_table;

  int core_start;
  int core_end;

  void (*coarse_hopp)(vector_PRECISION eta, vector_PRECISION phi, OPERATOR_TYPE_PRECISION *D, level_struct *l);
  if(sign == +1)
    coarse_hopp = coarse_hopp_PRECISION_vectorized;
  else
    coarse_hopp = coarse_n_hopp_PRECISION_vectorized;


  if ( l->num_processes > 1 && op->c.comm ) {
    set_boundary_PRECISION( out, 0, l, threading );

    if ( amount == _EVEN_SITES ) {
      minus_dir_param = _ODD_SITES;
      plus_dir_param = _EVEN_SITES;
    } else if ( amount == _ODD_SITES ) {
      minus_dir_param = _EVEN_SITES;
      plus_dir_param = _ODD_SITES;
    }

    START_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      // send in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
    END_MASTER(threading)

    if ( amount == _EVEN_SITES ) {
      start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
    } else if ( amount == _ODD_SITES ) {
      start = 0; num_lattice_sites = op->num_even_sites;
    }
    compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

    // prepare for sending to fw: compute hopping terms into forward boundary buffer
    for ( i=core_start; i<core_end; i++ ) {
      for(int mu=0; mu<4; mu++) {
        if(neighbor_fw[5*i+1+mu] < l->num_inner_lattice_sites)
          continue;
        out_pt = out + num_site_var*neighbor_fw[5*i+1+mu];
        in_pt = in + num_site_var*neighbor_fw[5*i];
        D_vectorized = op->D_transformed_vectorized + 4*link_offset*neighbor_fw[5*i] + mu*link_offset;
        coarse_hopp( out_pt, in_pt, D_vectorized, l );
      }
    }

    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      // send in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
    END_LOCKED_MASTER(threading)

  }
  else
    SYNC_CORES(threading)


  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  for ( i=core_start; i<core_end; i++ ) {
    out_pt = out + num_site_var*neighbor_fw[5*i];

    // U_mu^dagger coupling
    for(int mu=0; mu<4; mu++) {
      // terms coming from backward boundary buffer are done by the ghost_wait_PRECISION call below
      if(neighbor_bw[5*i+1+mu] >= l->num_inner_lattice_sites)
        continue;
      D_vectorized = op->D_transformed_vectorized + 4*link_offset*neighbor_bw[5*i+1+mu] + mu*link_offset;
      in_pt = in + num_site_var*neighbor_bw[5*i+1+mu];
      coarse_hopp( out_pt, in_pt, D_vectorized, l );
    }

    // compute U_mu couplings
    for(int mu=0; mu<4; mu++) {
      D_vectorized = op->D_vectorized + 4*link_offset*neighbor_fw[5*i] + mu*link_offset;
      in_pt = in + num_site_var*neighbor_fw[5*i+1+mu];
      coarse_hopp( out_pt, in_pt, D_vectorized, l );
    }
  }

  // wait for terms from bw and add them
  if ( l->num_processes > 1 && op->c.comm ) {

    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    END_LOCKED_MASTER(threading)

  }
  else
    SYNC_CORES(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_n_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading ) {

#ifdef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;

  OPERATOR_TYPE_PRECISION *D_vectorized;
  int column_offset = 2*SIMD_LENGTH_PRECISION*((l->num_parent_eig_vect+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int vectorized_link_offset = 2*2*l->num_parent_eig_vect*column_offset;

  int core_start;
  int core_end;

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );

  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)

  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // D is applied in an input-centric way
  // this makes threading a bit ugly, is there a better way?
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 0*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 1*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 2*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 3*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_vectorized + 4*vectorized_link_offset*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

  int fgmres_iters=0;

  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( p->b, p->x, op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );

  int start, end;
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

#ifdef BLOCK_JACOBI
  if ( l->level==0 && l->p_PRECISION.block_jacobi_PRECISION.local_p.polyprec_PRECISION.update_lejas == 1 ) {
    // re-construct Lejas
    local_re_construct_lejas_PRECISION( l, threading );
  }
#endif

#ifdef BLOCK_JACOBI
  // if Block Jacobi is enabled, solve the problem : M^{-1}Ax = M^{-1}b
  if ( p->block_jacobi_PRECISION.BJ_usable == 1 ) {
    // create a backup of b
    vector_PRECISION_copy( p->block_jacobi_PRECISION.b_backup, p->b, start, end, l );
    block_jacobi_apply_PRECISION( p->b, p->block_jacobi_PRECISION.b_backup, p, l, threading );
  }
#endif

#ifdef POLYPREC
  if ( l->level==0 && l->p_PRECISION.polyprec_PRECISION.update_lejas == 1 ) {
    // re-construct Lejas
    re_construct_lejas_PRECISION( l, threading );
  }
#endif

#ifdef POLYPREC
  START_MASTER(threading)
  p->preconditioner = p->polyprec_PRECISION.preconditioner;
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
#endif

#ifdef GCRODR
  fgmres_iters = flgcrodr_PRECISION( p, l, threading );
#else
  fgmres_iters = fgmres_PRECISION( p, l, threading );
#endif

  START_MASTER(threading)
  g.avg_b1 += fgmres_iters;
  g.avg_b2 += 1;
  g.avg_crst = g.avg_b1/g.avg_b2;
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

#ifdef BLOCK_JACOBI
  // restore the rhs
  if ( p->block_jacobi_PRECISION.BJ_usable == 1 ) {
    vector_PRECISION_copy( p->b, p->block_jacobi_PRECISION.b_backup, start, end, l );
  }
#endif

  SYNC_MASTER_TO_ALL(threading)

  START_MASTER(threading)
  if (g.my_rank==0) printf("coarsest gmres iters = %d\n", fgmres_iters);
  END_MASTER(threading)

  // even to odd
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( p->b, p->x, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}


void coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start, &end, l, threading);

  vector_PRECISION *tmp = op->buffer;

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );

  SYNC_CORES(threading)
  vector_PRECISION_define( tmp[0], 0, start, end, l );

#ifdef PERS_COMMS
  g.use_pers_comms2 = 1;
#endif
  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  //START_MASTER(threading)
  //printf0("print 1 ...\n");
  //END_MASTER(threading)
  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
#ifdef PERS_COMMS
  g.use_pers_comms2 = 0;
#endif

  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( tmp[1], tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );

#ifdef PERS_COMMS
  g.use_pers_comms2 = 1;
#endif
  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( out, tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
#ifdef PERS_COMMS
  g.use_pers_comms2 = 0;
#endif

}


void g5D_coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start_even, end_even, start_odd, end_odd;
  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );
  
  vector_PRECISION tmp = op->buffer[0];
  
  SYNC_CORES(threading)
  vector_PRECISION_define( tmp, 0, start_even, end_even, l );
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_gamma5_PRECISION( p->b, p->b, start_odd, end_odd, l );
  SYNC_CORES(threading)
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( p->b, p->b, start_odd, end_odd, l );
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( tmp, p->x, op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  coarse_gamma5_PRECISION( tmp, tmp, start_even, end_even, l );
  SYNC_CORES(threading)
  vector_PRECISION_plus( p->b, p->b, tmp, start_even, end_even, l );

#ifdef GCRODR
  flgcrodr_PRECISION( p, l, threading );
#else
  fgmres_PRECISION( p, l, threading );
#endif
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( p->b, p->b, start_odd, end_odd, l );
  SYNC_CORES(threading)
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  SYNC_CORES(threading)
  
  // even to odd
  PROF_PRECISION_START( _NC, threading );
  vector_PRECISION_define( tmp, 0, start_odd, end_odd, l );
  SYNC_CORES(threading)
  coarse_n_hopping_term_PRECISION( tmp, p->x, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->b, tmp, op, l, threading );
  vector_PRECISION_plus( p->x, p->x, p->b, start_odd, end_odd, l );
  
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}


void g5D_coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

  int start_even, end_even, start_odd, end_odd;
  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );

  vector_PRECISION *tmp = op->buffer;

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );

  vector_PRECISION_define( tmp[0], 0, start_odd, end_odd, l );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( tmp[1], tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( out, tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  coarse_gamma5_PRECISION( out, out, start_even, end_even, l );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

}


void coarse_odd_even_PRECISION_test( vector_PRECISION out, vector_PRECISION in, level_struct *l, struct Thread *threading ) {
  
  if ( g.odd_even ) {
    vector_PRECISION buf1 = NULL, buf2 = NULL;
    
    PUBLIC_MALLOC( buf1, complex_PRECISION, 2*l->vector_size );
    buf2 = buf1 + l->vector_size;

    START_LOCKED_MASTER(threading)
    // transformation part
    vector_PRECISION_copy( buf1, in, 0, l->inner_vector_size, l );
    // even to odd
    vector_PRECISION_define( out, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)

    coarse_hopping_term_PRECISION( out, buf1, &(l->oe_op_PRECISION), _ODD_SITES, l, threading );
    coarse_diag_oo_inv_PRECISION( buf2, out, &(l->oe_op_PRECISION), l, threading );

    START_LOCKED_MASTER(threading)
    vector_PRECISION_plus( buf1, buf1, buf2, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    
    // block diagonal part
    if ( g.method == 6 ) {
      g5D_coarse_apply_schur_complement_PRECISION( out, buf1, &(l->oe_op_PRECISION), l, threading );
    } else {
      coarse_apply_schur_complement_PRECISION( out, buf1, &(l->oe_op_PRECISION), l, threading );
    }
    
    coarse_diag_oo_PRECISION( out, buf1, &(l->oe_op_PRECISION), l, threading );
    
    // back transformation part
    coarse_diag_oo_inv_PRECISION( buf2, out, &(l->oe_op_PRECISION), l, threading );
    
    if ( g.method == 6 ) {
      START_LOCKED_MASTER(threading)
      coarse_gamma5_PRECISION( out, out, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
      vector_PRECISION_define( buf1, 0, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      coarse_hopping_term_PRECISION( buf1, buf2, &(l->oe_op_PRECISION), _EVEN_SITES, l, no_threading );
      coarse_gamma5_PRECISION( buf1, buf1, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      vector_PRECISION_plus( out, out, buf1, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      END_LOCKED_MASTER(threading)
    } else {
      coarse_hopping_term_PRECISION( out, buf2, &(l->oe_op_PRECISION), _EVEN_SITES, l, threading );
    }

    PUBLIC_FREE( buf1, complex_PRECISION, 2*l->vector_size );
  }
}

void coarse_apply_oddeven_operator_PRECISION(vector_PRECISION out,
                                                vector_PRECISION in,
                                                operator_PRECISION_struct *op,
                                                level_struct *l,
                                                struct Thread *threading){
    /* applies coarsest operator using different function calls (4 in total)
D = /D_ee D_eo\
    \D_oe D_oo/ 

    out_e = D_ee in_e + D_eo in_o
    out_o = D_oo in_o + D_oe in_e

*/

  //sets output to 0s
  vector_PRECISION_define(out, 0, l->p_PRECISION.v_start, 2* l->p_PRECISION.v_end, l);

  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start, &end, l, threading);

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading ); // out_e = D_ee in_e
  PROF_PRECISION_STOP( _SC, 0, threading );

  SYNC_CORES(threading)

  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( out, in, op, _EVEN_SITES, l, threading ); //out += D_eo in_o	will touch only the part out_e
  PROF_PRECISION_STOP( _NC, 1, threading );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_PRECISION( out, in, op, l, threading ); // out = D_oo in_o will touch only out_o
  PROF_PRECISION_STOP( _SC, 0, threading );

  SYNC_CORES(threading)

  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( out, in, op, _ODD_SITES, l, threading ); //out += D_oe in_e will touch only out_o
  PROF_PRECISION_STOP( _NC, 1, threading );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
}
