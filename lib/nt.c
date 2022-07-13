/*

nt.c - number theoretic helper functions

Copyright (C) 2009, 2010, 2015, 2021, 2022 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include "cm-impl.h"

static void tree_gcd (mpz_t *gcd, mpz_srcptr n, mpz_t *m, int no_m);
static int miller_rabin (mpz_srcptr n);

/*****************************************************************************/

long int cm_nt_gcd (long int a, long int b)
   /* returns the positive gcd of a and b */

{
   long int r, a_local, b_local;

   if (a == 0)
      return b;
   else if (b == 0)
      return a;
   else
   {
      if (a < 0)
         a_local = -a;
      else
         a_local = a;
      if (b < 0)
         b_local = -b;
      else
         b_local = b;
      r = a_local % b_local;
      while (r > 0)
      {
         a_local = b_local;
         b_local = r;
         r = a_local % b_local;
      }
      return b_local;
   }
}

/*****************************************************************************/

int cm_nt_kronecker (int_cl_t a, int_cl_t b)
   /* Compute the Kronecker symbol (a / b) following Algorithm 1.4.12
      in Cohen93 (by the binary algorithm). */
{
   const int tab [8] = {0, 1, 0, -1, 0, -1, 0, 1};
   int k;
   int_cl_t r;

   /* step 1 */
   if (b == 0) {
      if (a == 1 || a == -1)
         return 1;
      else
         return 0;
   }

   /* step 2*/
   if (a % 2 == 0 && b % 2 == 0)
      return 0;

   while (b % 4 == 0)
      b >>= 2;
   if (b % 2 == 0) {
      b >>= 1;
      k = tab [cm_classgroup_mod (a, (uint_cl_t) 8)];
   }
   else
      k = 1;

   if (b < 0) {
      b = -b;
      if (a < 0)
         k = -k;
   }

   /* step 3 */
   a = cm_classgroup_mod (a, b);

   /* steps 4 to 6 */
   while (a != 0) {
      while (a % 4 == 0)
         a >>= 2;
      if (a % 2 == 0) {
         a >>= 1;
         k *= tab [cm_classgroup_mod (b, (uint_cl_t) 8)];
      }
      if (b > a) {
         r = b - a;
         if (   cm_classgroup_mod (a, (uint_cl_t) 4) == 3
             && cm_classgroup_mod (b, (uint_cl_t) 4) == 3)
            k = -k;
         b = a;
         a = r;
      }
      else
         a -= b;
   }

   if (b > 1)
      return 0;
   else
      return k;
}

/*****************************************************************************/

static int miller_rabin (mpz_srcptr n)
   /* Return whether the odd positive integer n is a strong pseudoprime
      to base 2. */

{
   int d, res;
   mpz_t nm1, b, e;

   mpz_init (nm1);
   mpz_init (b);
   mpz_init (e);

   /* Write n-1 = 2^d * e. */
   mpz_sub_ui (nm1, n, 1);
   d = mpz_scan1 (nm1, 0);
   mpz_tdiv_q_2exp (e, nm1, d);

   mpz_set_ui (b, 2);
   mpz_powm (b, b, e, n);
   if (!mpz_cmp_ui (b, 1) || !mpz_cmp (b, nm1))
      res = 1;
   else {
      res = 0;
      while (!res && d > 1) {
         mpz_powm_ui (b, b, 2, n);
         d--;
         if (!mpz_cmp (b, nm1))
            res = 1;
      }
   }

   mpz_clear (nm1);
   mpz_clear (b);
   mpz_clear (e);

   return res;
}

/*****************************************************************************/

int cm_nt_is_prime (mpz_srcptr n)

{
   /* According to Table 1 of
      https://doi.org/10.1090/S0025-5718-1989-0982368-4
      the probability that a Miller-Rabin test fails for a random number
      and a random base for numbers with 900 digits is less than 10^(-109).
      We hope that this still holds for the fixed base 2. */
   if (mpz_sizeinbase (n, 2) >= 3000)
      return miller_rabin (n);
   else
      return (mpz_probab_prime_p (n, 0) > 0);
}

/*****************************************************************************/

unsigned long int cm_nt_next_prime (const unsigned long int n)
   /* returns the prime following n */

{
   static bool init = true;
#ifdef WITH_MPI
   static unsigned long int P [664579];
      /* primes up to 10^7 */
#else
   static unsigned long int P [9592];
      /* primes up to 10^5 */
#endif
   const unsigned int size = sizeof (P) / sizeof (unsigned long int);

   if (init) {
      mpz_t p;
      unsigned int i;

      mpz_init_set_ui (p, 0ul);
      for (i = 0; i < size; i++) {
         mpz_nextprime (p, p);
         P [i] = mpz_get_ui (p);
      }

      mpz_clear (p);
      init = false;
   }

   if (n < P [size - 1]) {
      /* search in prime list; loop invariants:
         left <= right
         P [0], ..., P [left - 1] <= n
         P [right], ..., P [size - 1] > n      */
      unsigned int left = 0, right = size - 1, middle;
      while (left != right) {
         middle = (left + right) / 2; /* < right */
         if (P [middle] > n)
            right = middle; /* becomes smaller */
         else
            left = middle + 1; /* becomes bigger, remains <= right */
      }
      return P [left];
   }
   else {
      printf ("***** Error: cm_nt_next_prime called with an argument\n"
         "that is too large for the precomputed list.\n");
      exit (1);
   }
}

/*****************************************************************************/

void cm_nt_factor (uint_cl_t d, uint_cl_t *factors, unsigned int *exponents)
   /* Factor d by trial division. The prime factors are stored in factors,
      their multiplicities in exponents, which must provide sufficiently
      much space. The list of prime factors is terminated by 0, so that
      17 entries suffice for the 64 bit type. */

{
   uint_cl_t p, p2;
   int i;

   i = 0;
   p = 2;
   p2 = 4;
   do {
      if (d % p == 0) {
         factors [i] = p;
         exponents [i] = 1;
         d /= p;
         while (d % p == 0) {
            exponents [i]++;
            d /= p;
         }
         i++;
      }
      /* We may wish to use cm_nt_next_prime, but its implementation
         is slow beyond the precomputed table. */
      if (p == 2) {
         p = 3;
         p2 = 9;
      }
      else {
         p2 += 4*(p+1);
         p += 2;
      }
   } while (p2 <= d);

   if (d != 1) {
     /* There is a prime factor left. */
     factors [i] = d;
     exponents [i] = 1;
     i++;
   }

   factors [i] = 0;
}

/*****************************************************************************/

uint_cl_t cm_nt_largest_factor (uint_cl_t n)
   /* Return the largest prime factor of n or 1 if n==1. */
{
   uint_cl_t factors [17];
   unsigned int exponents [17];
   int i;

   cm_nt_factor (n, factors, exponents);
   for (i = 0; factors [i] != 0; i++);

   if (i == 0)
      return 1;
   else
      return factors [i - 1];
}

/*****************************************************************************/

unsigned int cm_nt_mpz_tonelli_generator (mpz_ptr q, mpz_ptr z,
   mpz_srcptr p)
   /* For p an odd prime, compute and return e such that p-1 = 2^e * q with
      q odd. If e>1 (that is, p != 3 mod 4), also q and a generator z of
      the 2-Sylow subgroup of (Z/pZ)^* are returned. Otherwise q and z
      are not modified. */
{
   unsigned int e = mpz_scan1 (p, 1);

   if (e > 1) {
      mpz_tdiv_q_2exp (q, p, e);
      for (mpz_set_ui (z, 2ul); mpz_legendre (z, p) != -1;
         mpz_add_ui (z, z, 1ul));
      mpz_powm (z, z, q, p);
   }

   return e;
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_with_generator (mpz_ptr root, mpz_srcptr a,
   mpz_srcptr p, unsigned int e, mpz_srcptr q, mpz_srcptr z)
   /* Compute a square root of a modulo p by the Tonelli-Shanks algorithm,
      see Cohen93, Algorithm 1.5. */
{
   mpz_t a_local, y, x, b, tmp;
   unsigned long int r, m;

   mpz_init (a_local);
   mpz_tdiv_r (a_local, a, p);
   if (mpz_cmp_ui (a_local, 0ul) == 0) {
      mpz_set_ui (root, 0ul);
      mpz_clear (a_local);
      return;
   }

   mpz_init (y);
   mpz_init (x);
   mpz_init (b);
   mpz_init (tmp);

   if (e == 1) /* p=3 (mod 4) */ {
      mpz_add_ui (tmp, p, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 2ul);
      mpz_powm (x, a_local, tmp, p);
   }
   else {
      /* initialisation */
      mpz_set (y, z);
      r = e;
      mpz_sub_ui (tmp, q, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 1ul);
      mpz_powm (x, a_local, tmp, p);
      mpz_powm_ui (b, x, 2ul, p);
      mpz_mul (b, b, a_local);
      mpz_mod (b, b, p);
      mpz_mul (x, x, a_local);
      mpz_mod (x, x, p);
      while (mpz_cmp_ui (b, 1ul) != 0) {
         /* find exponent */
         m = 1;
         mpz_powm_ui (tmp, b, 2ul, p);
         while (mpz_cmp_ui (tmp, 1ul) != 0) {
            m++;
            mpz_powm_ui (tmp, tmp, 2ul, p);
         }
         if (m == r) {
            printf ("*** mpz_tonelli called with a = ");
            mpz_out_str (stdout, 10, a);
            printf (" and p = ");
            mpz_out_str (stdout, 10, p);
            printf (",\nbut a is not a square modulo p!\n");
            exit (1);
         }
         r = 1 << (r - m - 1);
         mpz_powm_ui (tmp, y, r, p);
         mpz_powm_ui (y, tmp, 2ul, p);
         r = m;
         mpz_mul (x, x, tmp);
         mpz_mod (x, x, p);
         mpz_mul (b, b, y);
         mpz_mod (b, b, p);
      }
   }

   mpz_set (root, x);

   mpz_clear (a_local);
   mpz_clear (y);
   mpz_clear (x);
   mpz_clear (b);
   mpz_clear (tmp);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli (mpz_ptr root, mpz_srcptr a, mpz_srcptr p)
{
   unsigned int e;
   mpz_t q, z;

   mpz_init (q);
   mpz_init (z);

   e = cm_nt_mpz_tonelli_generator (q, z, p);
   cm_nt_mpz_tonelli_with_generator (root, a, p, e, q, z);

   mpz_clear (q);
   mpz_clear (z);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_si_with_generator (mpz_ptr root,
   const long int a, mpz_srcptr p, unsigned int e, mpz_srcptr q,
   mpz_srcptr z)
{
   mpz_t tmp_a;

   mpz_init_set_si (tmp_a, a);
   cm_nt_mpz_tonelli_with_generator (root, tmp_a, p, e, q, z);
   mpz_clear (tmp_a);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_si (mpz_ptr root, const long int a, mpz_srcptr p)
   /* computes a square root of a modulo p */

{
   mpz_t tmp_a;

   mpz_init_set_si (tmp_a, a);
   cm_nt_mpz_tonelli (root, tmp_a, p);
   mpz_clear (tmp_a);
}

/*****************************************************************************/

static void tree_gcd (mpz_t *gcd, mpz_srcptr n, mpz_t *m, int no_m)
   /* This is the workhorse behind the cm_nt_mpz_tree_gcd function; it
      assumes that the total size of the no_m entries of m is less than
      that of n and directly computes a tree. */
{
   mpz_t **tree;
   int *width;
   int levels;
   int i, j;

   /* Compute the height of the subproduct tree of the m in levels. */
   for (i = no_m, levels = 1; i > 1; i = (i+1) / 2, levels++);

   /* Compute bottom-up a subproduct tree with m on the leaves. */
   tree = (mpz_t **) malloc (levels * sizeof (mpz_t *));
   width = (int *) malloc (levels * sizeof (int));
   width [0] = no_m;
   tree [0] = (mpz_t *) malloc (no_m * sizeof (mpz_t));
   for (j = 0; j < no_m; j++)
      mpz_init_set (tree [0][j], m [j]);
   for (i = 1; i < levels; i++) {
      width [i] = (width [i-1] + 1) / 2;
      tree [i] = (mpz_t *) malloc (width [i] * sizeof (mpz_t));
      for (j = 0; j < width [i-1] / 2; j++) {
         mpz_init (tree [i][j]);
         mpz_mul (tree [i][j], tree [i-1][2*j], tree [i-1][2*j+1]);
      }
      if (width [i-1] % 2 != 0) {
         mpz_init (tree [i][j]);
         mpz_set (tree [i][j], tree [i-1][2*j]);
      }
   }

   /* Replace the tree tops by n modulo the entry. */
   for (j = 0; j < width [levels-1]; j++)
      mpz_mod (tree [levels-1][j], n, tree [levels-1][j]);

   /* Replace top-down the tree entries by n modulo the entry. */
   for (i = levels - 2; i >= 0; i--) {
      for (j = 0; j < (width [i] / 2) * 2; j++)
         mpz_mod (tree [i][j], tree [i+1][j/2], tree [i][j]);
      if (width [i] % 2 != 0)
         mpz_set (tree [i][j], tree [i+1][j/2]);
   }

   /* Compute the gcd of n mod m [j] and m [j]. */
   for (j = 0; j < no_m; j++)
      mpz_gcd (gcd [j], tree [0][j], m [j]);

   /* Clear the tree. */
   for (i = 0; i < levels; i++) {
      for (j = 0; j < width [i]; j++)
         mpz_clear (tree [i][j]);
      free (tree [i]);
   }
   free (tree);
   free (width);
}

/*****************************************************************************/

void cm_nt_mpz_tree_gcd (mpz_t *gcd, mpz_srcptr n, mpz_t *m, int no_m)
   /* Given a (large) positive integer n and an array of (smaller) no_m >= 1
      positive integers in m, compute the gcd of n with all the m and return
      them in gcd, which needs to provide sufficient space and initialised
      entries. */
{
   const unsigned long int size_n = mpz_sizeinbase (n, 2);
   unsigned long int size_m;
   int offset, no_batch;

   for (offset = 0; offset < no_m; offset += no_batch) {
      /* Treat the entries m [offset], ..., m [offset + no_batch - 1],
         where no_batch is chosen such that their product is at most
         size_n / 2; this keeps the subproduct tree of the m, which
         can require a considerable amount of memory, of a size such
         that reduction of n modulo the product at the root still
         requires some work. */
      for (size_m = mpz_sizeinbase (m [offset], 2), no_batch = 1;
           offset + no_batch < no_m
              && size_m + mpz_sizeinbase (m [offset + no_batch], 2)
                 < size_n / 2;
           size_m += mpz_sizeinbase (m [offset + no_batch], 2), no_batch++);
      tree_gcd (gcd + offset, n, m + offset, no_batch);
   }
}

/*****************************************************************************/
/*****************************************************************************/

bool cm_nt_fget_z (mpz_t out, ftype in)
   /* Tries to round the real value in to an integer value out. The return
      value reflects the success of the operation. */

{
   ftype rounded, diff;
   bool ok;
   mp_exp_t expo;

   finit (rounded, fget_prec (in));
   finit (diff, fget_prec (in));

   fround (rounded, in);
   fsub (diff, in, rounded);
   if (fsgn (diff) == 0 || (- fget_exp (diff)) >= 10) {
      expo = fget_z_exp (out, rounded);
      if (expo > 0)
         ok = false;
      else if (expo < 0) {
         mpz_tdiv_q_2exp (out, out, (unsigned long int) (-expo));
         ok = true;
      }
      else
         ok = true;
   }
   else
      ok = false;

   if (!ok) {
      printf ("***** Error in rounding:\n");
      fout_str (stdout, 10, 0ul, in);
      printf ("\n");
      fout_str (stdout, 10, 0ul, rounded);
      printf ("\n");
   }

   fclear (rounded);
   fclear (diff);

   return ok;
}

/*****************************************************************************/

bool cm_nt_cget_zz (mpz_ptr out1, mpz_ptr out2, ctype in, ctype omega)
   /* Tries to round the complex value to an imaginary-quadratic integer
      out1+out2*omega, where omega is the second element of an integral
      basis (or more generally, a non-real complex number). The return
      value reflects the success of the operation. */
{
   ftype tmp;
   bool ok;

   finit (tmp, cget_prec (in));

   fdiv (tmp, cimagref (in), cimagref (omega));
   ok = cm_nt_fget_z (out2, tmp);

   if (ok) {
      fmul_z (tmp, crealref (omega), out2);
      fsub (tmp, crealref (in), tmp);
      ok = cm_nt_fget_z (out1, tmp);
   }

   fclear (tmp);

   return ok;
}

/*****************************************************************************/
/*****************************************************************************/

void mpz_binary_prod(mpz_ptr ans,mpz_t* vec,unsigned long int size){
// set ans=product(i=0,size-1,vec[i])
    if(size==0){
        mpz_set_ui(ans,1);// empty product
        return;
    }
    
    if(size==1){
        mpz_set(ans,vec[0]);
        return;
    }
    
    unsigned long int mid=size>>1;
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_binary_prod(tmp,vec,mid);
    mpz_binary_prod(ans,vec+mid,size-mid);
    mpz_mul(ans,ans,tmp);
    mpz_clear(tmp);
}

/*****************************************************************************/

void mpz_binary_prod64(mpz_ptr ans,unsigned long int* vec,unsigned long int size){
// set ans=product(i=0,size-1,vec[i])
    if(size==0){
        mpz_set_ui(ans,1);// empty product
        return;
    }
    
    if(size==1){
        mpz_set_ui(ans,vec[0]);
        return;
    }
    
    unsigned long int mid=size>>1;
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_binary_prod64(tmp,vec,mid);
    mpz_binary_prod64(ans,vec+mid,size-mid);
    mpz_mul(ans,ans,tmp);
    mpz_clear(tmp);
}

/*****************************************************************************/

void cm_prime_product_basecase (mpz_ptr prim, unsigned long int a,  unsigned long int b, unsigned int* plist,unsigned int plist_size){
// set prim=prod(primes in [a,b) interval), not including b. 
// plist contains all primes up to sqrt(b), possibly more.
    if(a>=b || b<=2){
        mpz_set_ui(prim,1);
        return;
    }
    
    unsigned int j,*T;
    unsigned long int i,interval=sqrt(b)+1;
    unsigned long int num_primes=0;
    unsigned long int *vec;
    bool mult2=(a<=2&&2<b);
    
    if(a%2==0)a++;
    if(a<3)   a=3;
    interval=((interval+32)>>5)<<5;
    unsigned long int Tsize=interval>>5;
    T=(unsigned int*)malloc(Tsize*sizeof(unsigned int));
    
    unsigned long int maxsize=14.0+(double)1.25506*b/log(b)-(double)a/log(a);// x/log(x)<pi(x)<1.25506*x/log(x) for x>=17
    vec=(unsigned long int*)malloc(maxsize*sizeof(unsigned long int));
    
    for(;a<=b;a+=2*interval){ 
         for(i=0;i<Tsize;i++)T[i]=0;
         for(i=1;i<plist_size;i++){
             unsigned int p=plist[i];
             if(p*p>b)break;
             unsigned int res=a%p;
             if(res&1) res=(p-res)>>1;
             else      res=(2*p-res)>>1;
             if(res==p)res=0;
             if(a+2*res==p)res+=p;
             
             for(j=res;j<interval;j+=p)T[j>>5]|=(1u<<(j&31));
         }
         for(i=0;i<Tsize;i++){
             unsigned int tmp=T[i];
             for(j=0;j<32;j++,tmp>>=1)
                 if((tmp&1)==0){
                    unsigned long int q=a+64*i+2*j;
                    if(q<b)vec[num_primes++]=q;                     
                 }
         }
    }
    free(T);
    
    mpz_binary_prod64(prim,vec,num_primes);
    if(mult2)mpz_mul_ui(prim,prim,2);
    
    free(vec);
}

/*****************************************************************************/

void cm_prime_product_memoryefficient (mpz_ptr prim, unsigned long int a,  unsigned long int b){
//  set prim=prod(primes in [a,b) interval), not including b. 
    
    if(a>=b || b<=2) {
        mpz_set_ui(prim,1);
        return;
    }
    if (a<2) a=2;
    
    bool *T;
    unsigned int* plist;
    unsigned int i,j,plist_size=0;
    unsigned int sq=(unsigned int)sqrt(b)+2;
    
    T=(bool*)malloc((sq+1)*sizeof(bool));
    for(i=0;i<=sq;i++)T[i]=(i>1);
    for(i=2;i*i<=sq;i++)
        if(T[i]){
           for(j=i*i;j<=sq;j+=i)T[j]=false;
        }
    
    for(i=2;i<=sq;i++)if(T[i])plist_size++;
    plist=(unsigned int*)malloc(plist_size*sizeof(unsigned int));
    plist_size=0;
    for(i=2;i<=sq;i++)if(T[i])plist[plist_size++]=i;
    free(T);

#define maxdepth 4 // larger value decreases the memory but will increase the running time,
                   // we do the product in 2^maxdepth chunks.
    
    unsigned long int step=((b-a)>>maxdepth)+1;
        
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(prim,1);
    for(i=0;i<(1<<maxdepth);i++){
        unsigned long int a2=(a+step<=b?a+step:b);
        cm_prime_product_basecase(tmp,a,a2,plist,plist_size);
        mpz_mul(prim,prim,tmp);
        a+=step;
    }
    mpz_clear(tmp);      
    free(plist);
}

/*****************************************************************************/
/*****************************************************************************/