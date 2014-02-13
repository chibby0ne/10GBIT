#ifndef GAUSSIAN_RANDOM_GENERATOR_H
#define GAUSSIAN_RANDOM_GENERATOR_H

#include <math.h>
#ifdef __SSE2__
#include <emmintrin.h>
#endif

#include <time.h>
using namespace std;

/// Random
/** Gaussian_Random_Generator
 * used to generate either equaly or normally distributed numbers. 
 * uses simd if the -msse2 (SSE2) flag is set.
 */
class Gaussian_Random_Generator
{
	public:
		/// Standard constructor, initializes with a default seed
		Gaussian_Random_Generator() { if (!initialized) reset(4357U); }
		/// Initializes random numbers wirt a given seed
		Gaussian_Random_Generator(unsigned int seed) { reset(seed); }
		/// resets the random numbers with a new seed based on current time
		void randomize() { reset(hash(time(0), clock())); }
		/// resets the random numbers with the last seed
		void reset() { initialize(last_seed); reload(); initialized = true; }
		/// resets the random numbers with a given seed
		void reset(unsigned int seed) { last_seed = seed; reset(); }

		unsigned int random_int() {
			if (left == 0) reload();
			--left;

			unsigned int s1;
			s1 = *pNext++;
			s1 ^= (s1 >> 11);
			s1 ^= (s1 <<  7) & 0x9d2c5680U;
			s1 ^= (s1 << 15) & 0xefc60000U;
			return (s1 ^(s1 >> 18));
		}

		/**random_01
		 * returns a equally distributed random number in range of (0,1)
		 */
		double random_01() { return (random_int() + 0.5) * (1.0 / 4294967296.0); }
		/**random_01_lclosed
		 * returns a equally distributed random number in range of [0,1)
		 */
		double random_01_lclosed() { return random_int() * (1.0 / 4294967296.0); }
		/**random_01_closed
		 * returns a equally distributed random number in range of [0,1]
		 */
		double random_01_closed() { return random_int() * (1.0 / 4294967295.0); }
		double random53_01_lclosed() {
			return ((random_int() >> 5) * 67108864.0 + (random_int() >> 6))
				* (1.0 / 9007199254740992.0); // by Isaku Wada
		}


#ifdef __SSE2__
		/**
		 * RandomSIMD
		 * param: double rueckgabe[2], double rand1, double rand2
		 * turns the equally distributed rand1 and rand2 into normally distributed numbers which are stored in the given array
		 * uses simd (single instruction multiple data)
		 */
		void RandomSIMD(double rueckgabe[2], double rand1, double rand2);
#endif

		/**RandomSingle
		 * param: double rand
		 * turns the equally distributed rand into normally distibuted number and returns it
		 */
		double RandomSingle(double rand);

		/**get_random
		 * param: double output[2]
		 * generates two normally distributed random numbers and stors them in output. 
		 * uses smid if SSE2 flag is set
		 */
		void get_random(double output[2]);

	private:
		static bool initialized;
		unsigned int last_seed;
		static unsigned int state[624];
		static unsigned int *pNext;
		static int left;

#ifdef __SSE2__
		void beasleyBoth(double rueckgabe[2], __m128d x);
		void moroBoth(double rueckgabe[2], double d1, double d2, double x1, double x2);
#endif


		void initialize(unsigned int seed) {
			register unsigned int *s = state;
			register unsigned int *r = state;
			register int i = 1;
			*s++ = seed & 0xffffffffU;
			for (; i < 624; ++i) {
				*s++ = (1812433253U * (*r ^(*r >> 30)) + i) & 0xffffffffU;
				r++;
			}	 
		}

		void reload() {
			register unsigned int *p = state;
			register int i;
			for (i = 624 - 397; i--; ++p)
				*p = twist(p[397], p[0], p[1]);
			for (i = 397; --i; ++p)
				*p = twist(p[397-624], p[0], p[1]);
			*p = twist(p[397-624], p[0], state[0]);

			left = 624, pNext = state;
		}
		unsigned int hiBit(const unsigned int& u) const { return u & 0x80000000U; }
		unsigned int loBit(const unsigned int& u) const { return u & 0x00000001U; }
		unsigned int loBits(const unsigned int& u) const { return u & 0x7fffffffU; }
		unsigned int mixBits(const unsigned int& u, const unsigned int& v) const
		{ return hiBit(u) | loBits(v); }

		/*
		 * ----------------------------------------------------------------------
		 * --- ediap - 2007/01/17 ---
		 * ----------------------------------------------------------------------
		 * Wagners's implementation of the twist() function was as follows:
		 *  { return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfU); }
		 * However, this code caused a warning/error under MSVC++, because
		 * unsigned value loBit(s1) is being negated with `-' (c.f. bug report
		 * [1635988]). I changed this to the same implementation as is used in
		 * original C sources of Mersenne Twister RNG:
		 *  #define MATRIX_A 0x9908b0dfUL
		 *  #define UMASK 0x80000000UL
		 *  #define LMASK 0x7fffffffUL
		 *  #define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
		 *  #define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))
		 * ----------------------------------------------------------------------
		 */
		unsigned int twist(const unsigned int& m, const unsigned int& s0,
				const unsigned int& s1) const
		{ return m ^(mixBits(s0, s1) >> 1) ^(loBit(s1) ? 0x9908b0dfU : 0U); }
		unsigned int hash(time_t t, clock_t c);

#ifdef __SSE2__
		__m128d mul(__m128d v1, __m128d v2) {
			return _mm_mul_pd(v1, v2);
		}

		__m128d mul(__m128d v, double d) {
			__m128d s = _mm_set1_pd(d);
			return _mm_mul_pd(v, s);
		}

		__m128d add(__m128d v, double d) {
			__m128d s = _mm_set1_pd(d);
			return _mm_add_pd(v, s);
		}

		__m128d add(__m128d v1, __m128d v2) {
			return _mm_add_pd(v1, v2);
		}

		__m128d sub(__m128d v, double d) {
			__m128d s = _mm_set1_pd(d);
			return _mm_sub_pd(v, s);
		}

		__m128d sub(__m128d v1, __m128d v2) {
			return _mm_sub_pd(v1, v2);
		}

		__m128d div(__m128d v, double d) {
			__m128d s = _mm_set1_pd(d);
			return _mm_div_pd(v, s);
		}

		__m128d div(__m128d v1, __m128d v2) {
			return _mm_div_pd(v1, v2);
		}
#endif
};
#endif
