/* PageRank 
   The genetic.dat dataset comes from:
   http://www.cs.toronto.edu/~tsap/experiments/datasets/
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* allocate one object of given type */
#define NEW(type) ((type*)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define NEW_A(num,type) ((type*)calloc((size_t)(num),(size_t)sizeof(type)))

typedef unsigned int u_int;

/* vector */
typedef struct
{
  u_int dim;
  double *e;
} VEC;

/* row of a sparse matrix */
typedef struct
{
  u_int  nnz;  /* # of non-zero (nz) value on this row */
  u_int  *col; /* column identifier for each nz value */
  double *val; /* value for each nz value */
} SROW;

/* sparse matrix */
typedef struct
{
  u_int m, n;
  SROW *row;
} SMAT;

/* v_get -- gets a VEC of dimension 'dim'
   Precondition: size >= 0
   Postcondition: initialized to zero */
VEC *v_get( u_int size )
{
  VEC *v;
  
  if( (v = NEW(VEC)) == (VEC *)NULL )
  {
    fprintf( stderr, "v_get memory error" );
    exit( -1 );
  }
  
  v->dim = size;
  if( (v->e = NEW_A(size,double)) == (double *)NULL )
  {
    free( v );
    fprintf( stderr, "v_get memory error" );
    exit( -1 );
  }
  
  return v;
}

/* v_free -- returns VEC & associated memory back to memory heap */
int v_free( VEC *v )
{
  if( v == (VEC *)NULL )
    return -1;
  
  if( v->e == (double *)NULL ) 
  {
    free( v );
  }
  else
  {
    free( v->e );
    free( v );
  }
  
  return 0;
}

/* sm_get -- gets an mxn sparse matrix by dynamic memory allocation 
   Precondition: m>=0 && n>=0
   Postcondition: each row is empty*/
SMAT *sm_get( u_int m, u_int n )
{
  SMAT *M;
  u_int i;
  
  if( (M = NEW(SMAT)) == (SMAT *)NULL )
  {
    fprintf( stderr, "sm_get memory error" );
    exit( -1 );
  }
  
  M->m = m ; M->n = n;

  if( (M->row = NEW_A(m,SROW)) == (SROW *)NULL )
  {
    free( M );
    fprintf( stderr, "sm_get memory error" );
    exit( -1 );
  }
  
  for( i = 0 ; i < m ; i++ )
  {
    (M->row[i]).nnz = 0;
    (M->row[i]).col = (u_int *) NULL;
    (M->row[i]).val = (double *) NULL;
  }

  return M;
}

/* sm_free -- returns SMAT & associated memory back to memory heap */
int sm_free( SMAT *M )
{
  u_int i;
  SROW *ri;

  if( M == (SMAT *)NULL ) return -1;
  
  if( M->row == (SROW *)NULL ) 
  {
    free( M );
  }
  else
  {
    for( i = 0 ; i < M->m ; i++ )
    {
      ri = &(M->row[i]);
      if( ri->nnz > 0 )
      {
        free( ri->col );
        free( ri->val );
      }
    }
    free( M->row );
    free( M );
  }
  
  return 0;
}

/* sm_input -- file input of sparse matrix 
   Precondition: will only work with a binary matrix. */
SMAT *sm_input( FILE *fp )
{
  SMAT *M;
  u_int *col; /* temp array to store the nz col of the current row */
  u_int r;    /* index of current row */
  int c;      /* index of current column */
  SROW *ri;   /* pointer to the current row in M */
  u_int m,n,i,j,k;
  
  /* get dimension */
  if( fscanf( fp, " SparseMatrix: %u by %u", &m, &n ) < 2 )
  {
    fprintf( stderr, "sm_input error reading dimensions" );
    exit( -1 );
  }
  
  if( (col = NEW_A(n,u_int)) == (u_int *)NULL )
  {
    fprintf( stderr, "sm_input memory error" );
    exit( -1 );
  }

  M = sm_get( m, n );
  
  /* get entries */
  for( i=0 ; i<m ; i++ )
  {
    if( fscanf( fp, " row %u:", &r ) < 1 )
    {
      fprintf( stderr, "sm_input error reading line %u", i );
      exit( -1 );
    }
    ri = &(M->row[i]);
    j = 0;
    for( ; ; )
    {
      if( fscanf( fp, "%d", &c ) < 1 )
      {
        fprintf( stderr, "sm_input error reading line %u col x", i );
        exit( -1 );
      }
      if( c < 0 ) break;
      col[j] = c;
      j++;
    } /* j is the number of nz value in row i */

    if( ( (ri->col = NEW_A(j,u_int)) == (u_int *)NULL ) && ( j!=0 ) )
    {
      fprintf( stderr, "sm_input memory error" );
      exit( -1 );
    }
    if( ( (ri->val = NEW_A(j,double)) == (double *)NULL ) && ( j!=0 ) )
    {
      fprintf( stderr, "sm_input memory error" );
      exit( -1 );
    }

    ri->nnz = j;

    for( k = 0 ; k < j ; k++ )
    {
      ri->col[k] = col[k]; 
      ri->val[k] = 1.0;
    }
  }
  
  free( col );

  return M;
}

/* sm_output -- file output of sparse matrix 
   Postcondition: the result is not a valid entry for sm_input,
     since it also works for a non-binary matrix. */
void sm_output( FILE *fp, SMAT *M )
{
  u_int i,j;
  SROW *ri;

  fprintf( fp, "SparseMatrix: %d by %d\n", M->m, M->n );

  for( i = 0 ; i < M->m ; i++ )
  {
    fprintf( fp, "row %u: ", i ); 
    ri = &( M->row[i] );
    for( j = 0 ; j < ri->nnz ; j++ )
    {
      fprintf( fp, "%u:%1.5g ", ri->col[j], ri->val[j] );
    }
    fprintf( fp, "-1\n" );
  }
}

/* v_output -- file output of vector */
void v_output( FILE *fp, VEC *v )
{
  u_int i;

  fprintf( fp, "Vector: %d\n", v->dim );
  for( i = 0 ; i < v->dim ; i++ ) fprintf( fp, "%1.5g ", v->e[i] );
  putc( '\n', fp );
}


/* create_HM -- creation of HM by calcul of degree(i) for each M[i,j] nnz */
/* impact on SM */
SMAT* create_HM(SMAT* SM)
{
    u_int i,j;
    for(i=0;i<SM->m;++i)
    {
      for(j=0;j<SM->row[i].nnz;++j)
      {
        SM->row[i].val[j]=(SM->row[i].val[j])/(SM->row[i].nnz);
      }
    }

    return SM;
}

/* pr_first -- first implementation of the page rank algorithm */
VEC * pr_first(VEC * rk, SMAT * HM)
{

  VEC * rk_next = v_get(rk->dim);
  u_int i, j;

  for(i=0; i<rk_next->dim; ++i)
  {
    if(HM->row[i].nnz==0) //REPARTIR EQUITABLEMENT
    {
      for(j=0;j<rk_next->dim;++j)
      {
        rk_next->e[j]+=(rk->e[i])*(1.0/(rk_next->dim));
      }
    }
    else //ALGORITHME NORMAL
    {
      for(j=0; j<HM->row[i].nnz;++j)
      {
        rk_next->e[HM->row[i].col[j]]+=rk->e[i]*HM->row[i].val[j];
      }
    }
  }

  v_free(rk);
  return rk_next;
}

/* pr_second -- second implementation of the page rank algorithm with ergodicity*/
VEC * pr_second(VEC * rk, SMAT * HM, VEC* a, double alpha)
{

  double composanteSurfeur=0.0;
  VEC * rk_next = v_get(rk->dim);
  u_int i, j;

  for(i=0; i<rk->dim;++i)
  {
    composanteSurfeur+=alpha*(a->e[i])*(rk->e[i]);
  }

  composanteSurfeur+=(1-alpha);
  composanteSurfeur/=rk->dim;

  for(i=0; i<rk_next->dim; ++i)
  {    
    for(j=0; j<HM->row[i].nnz;++j)
    {
     rk_next->e[HM->row[i].col[j]]+=alpha*rk->e[i]*HM->row[i].val[j];
    }
  }

  for(i=0; i<rk_next->dim; ++i)
  {
    rk_next->e[i]+=composanteSurfeur;
  }

  v_free(rk);
  return rk_next;
}

int main()
{
  FILE *fp;
  SMAT *SM;
  SMAT *HM;

  fp = fopen( "exemple.dat", "r" );
  SM = sm_input( fp );
  fclose( fp );

  sm_output( stdout, SM );

  //Calcul de HM
  HM = create_HM(SM);

  sm_output( stdout, HM);

  //Exécution de l'implémentation 1 de PageRank
  VEC * rk = v_get(HM->n);
  u_int i;
  int n = 0;
  /*
  printf("\nPremière implémentation via matrice stochastique :\n");

  for(i=0;i<rk->dim;++i)
  {
    rk->e[i]=(1.0/(rk->dim));
  }

  printf("Itération %d : ",n++);
  v_output(stdout,rk);

  while(n<30)
  {
    rk=pr_first(rk,HM);
    printf("Itération %d : ",n++);
    v_output(stdout,rk);
  }
  */

  printf("\nSeconde implémentation avec ergodicité :\n");
  double alpha = 0.85;
  for(i=0;i<rk->dim;++i)
  {
    rk->e[i]=(1.0/(rk->dim));
  }

  VEC * a = v_get(HM->n);
  for(i=0;i<a->dim;++i)
  {
    if(HM->row[i].nnz==0)
    {
      a->e[i]=1;
    }
    else
    {
      a->e[i]=0;
    }
  }

  printf("Vecteur a (noeuds absorbants) : ");
  v_output(stdout, a);

  n=0;
  while(1)
  {
    rk=pr_second(rk,HM,a,alpha);
    printf("Itération %d : ",n++);
    v_output(stdout,rk);
  }

  v_free(rk);
  sm_free( SM );

  return 0;
}
