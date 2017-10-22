
/* compute number of polygons with integer sides and perimeter n
   for n=3 to 250

   Tested on Ubuntu Linux 16.04 LTS 64-bit.
   Requires installation of boost package, 
      sudo apt-get install libboost-all-dev
   Compile: 
      g++ -O3 -o mgons mgons.cpp
*/

#include <iostream>
#include <stdexcept>

#include <boost/multiprecision/cpp_int.hpp>
typedef boost::multiprecision::int256_t multiprecision;

/* Euclidean GCD algorithm */
unsigned compute_gcd(unsigned a, unsigned b) {
   for(;;) {
      unsigned z = b % a;
      if (z == 0) break;
      b = a;
      a = z;
   }
   return a;
}

/* Euler's totient: # of natural numbers less than and relatively prime to n */
unsigned totient(unsigned z) {
   unsigned res = 1;
   for (unsigned i=2; i<z; ++i) {
      if (compute_gcd(i,z) == 1) {
         ++res;
      }
   }
   return res;
}

/* binomial coefficent "bn choose bd" These can get big so use multiprecision integers */
multiprecision binom(unsigned bn, unsigned bd) {
   multiprecision res = 1;
   for (unsigned i=1; i<=bd; ++i) {
      res *= bn--;
      res /= i;
   }
   return res;
}

/* Compute the number of m-gons with integer sides and perimeter n, using formula from East-Niles paper */

multiprecision compute_mgon(unsigned m, unsigned n) {
   unsigned denom = 4* m; /* want to work with integers, so this is the worst-case LCD */
   multiprecision num = 0;
   unsigned gcd = compute_gcd(m, n);

   /* first term: summation over divisors of m and n */
   for (unsigned d=1; d<=gcd; ++d) {
      if (gcd % d != 0) continue;
      num += 2*totient(d)*binom(n / d - 1, m / d - 1); /* multiply by 2, it will become 1/2 when we divide by "denom" */
   }

   /* Second term, multiply by 2*m, which will become 1/2 when we divide by "denom" */
   num -= 2* m * binom(n/2, m - 1);

   /* now, separate cases for m even and odd */
   if (m & 1) {
      /* odd case */
      multiprecision msum = 0;
      for (unsigned i=1; i<=(n-1)/2; ++i) {
         if ((n+i)&1) continue; /* only compute when "perimiter" and "i" have the same parity */
         msum += binom((n-i)/2-1, m/2-1);
      }
      num += msum * m *2; /* multiply by m * 2, which will become 1/2 when we divide by "denom" */
   } else {
      /* even case */
      multiprecision msum = 0;

      /* if n is also even we need this term */
      if ((n&1)==0) {
         msum += binom(n/2-1, m/2-1);
      }
      /* then this double summation, whenever n-i-j is even */
      for (unsigned i=1; i<=(n-1)/2; ++i) {
         for (unsigned j=1; j<=(n-1)/2; ++j) {
            unsigned bn = n - i - j;
            if (bn&1) continue;
            msum += binom(bn/2-1, m/2-2);
         }
      }
      num += msum * m; /* multiply by m, this will become 1/4 when we divide by "denom" */
   }

   /* "sanity check:" Burnside's Lemma guarantees this evaluates to an integer, so throw an exception
      if this is not the case */
   if (num % denom != 0) throw std::runtime_error("bad_fraction");
   return num / denom;
}

int main(int argc, char *argv[]) {
   try {
      /* compute up to n=250, this is close to the limit of a 256-bit integer */
      for (unsigned i=3; i<=250; ++i) {
         /* add up all of the possible m-gons */
         multiprecision sum = 0;
         for (unsigned j=3; j<=i; ++j) {
            sum += compute_mgon(j,i);
         }
         /* compute a pseudo-log2 for the result, it's actually the ceiling of log2 */
         unsigned log2 = 0;
         multiprecision ml = 1;
         while(log2 < 256) {
            if (ml >= sum) break;
            ++log2;
            ml <<=1;
         }
         /* print a few comma separated lines, put quotes around big number so we have option to treat it as a string 
            when importing to spreadsheet */
         std::cout << "\"" << sum << "\"," << i << "," << log2 << std::endl;
      }
   } catch (std::exception &x) {
      std::cout << std::endl << x.what() << std::endl;
      return 1;
   }
   return 0;
}
