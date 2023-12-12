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

// vector storage for PRECISION precision
void vector_PRECISION_define( vector_PRECISION phi, complex_PRECISION value, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = value;
  } else {
    error0("Error in \"vector_PRECISION_define\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_define_random( vector_PRECISION phi, int start, int end, level_struct *l ) {

  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = (PRECISION)(((double)rand()/(double)RAND_MAX))-0.5 + ( (PRECISION)((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I;
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_define_random_rademacher( vector_PRECISION phi, int start, int end, level_struct *l ) {

  // TODO : use l->use_dilution to implement dilution or not
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ ) {

      //if ( g.use_dilution[l->depth]==0 ) {
      if ( l->use_dilution==0 ) {
        if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i]=  (double) (-1);
        else phi[i]= (PRECISION)(1);
      } else {
        if ( l->dil_spin==0 ) {
          if ( i%12 < 6 ) {
            if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i]=  (double) (-1);
            else phi[i]= (PRECISION)(1);
          }
          else {
            phi[i] = (PRECISION)0.0;
          }
        }
        else {
          if ( i%12 < 6 ) {
            phi[i] = (PRECISION)0.0;
          }
          else {
            if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i]=  (double) (-1);
            else phi[i]= (PRECISION)(1);
          }
        }
      }

    }
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );
  
}
