#include "vector_mod.h"
#include "mod_ops.h"
#include <vector>
#include <thread>
#include <cassert>
#include "num_threads.h"

#if defined(__GNUC__) && __GNUC__ < 10
namespace std
{
	constexpr std::size_t hardware_constructive_interference_size = 64u;
	constexpr std::size_t hardware_destructive_interference_size = 64u; //See False sharing
}
#endif

#ifndef NDEBUG
#define assert_(x) do {if (!(x)) {fprintf(stderr, "Error in %s (%s:%d)\n", #x, __FILE__, __LINE__); std::abort();}} while(false)
#else
#define assert_(x) 
#endif

constexpr std::size_t C = std::hardware_constructive_interference_size;
//const std::size_t T = std::thread::hardware_concurrency();

static std::vector<IntegerWord> get_powers_of_w(IntegerWord mod)
{
	std::vector<IntegerWord> result;
	std::size_t T = get_num_threads();
	result.reserve(C + T);
	result.emplace_back(1);
	for (std::size_t c = 1; c <= C; ++c)
		result.emplace_back(times_word(result.back(), mod));
	for (std::size_t t = 2; t <= T; ++t)
		result.emplace_back(mul_mod(result.back(), result[C], mod));
	return result;
}

static IntegerWord reduce_subvector(const IntegerWord* V, std::size_t n, IntegerWord mod,
	const IntegerWord* powers_of_w)
{
	assert(n <= C);
	IntegerWord result = 0;
	for (std::size_t c = 0; c < n; ++c)
		result = add_mod(result, mul_mod(V[c], powers_of_w[c], mod), mod);
	return result;
}

struct partial_result_t
{
	alignas(std::hardware_destructive_interference_size) IntegerWord value;
};

IntegerWord vector_mod(const IntegerWord* V, std::size_t N, IntegerWord mod)
{
	std::vector<IntegerWord> powers_of_w = get_powers_of_w(mod);
	std::size_t T = get_num_threads();
	std::vector<partial_result_t> partial_results(T);
	auto thread_proc = [&powers_of_w, &partial_results, V, N, mod, T](unsigned t)
	{
		IntegerWord result = 0;
		std::size_t subvector_count = ceil_div(N, C);
		std::size_t s_max = subvector_count - (subvector_count % T) + t;
		if (s_max >= subvector_count)
		{
			if (s_max < T)
				return;
			s_max -= T;
		}
		for (std::size_t s = s_max; true; s -= T)
		{
			std::size_t subvector_size = C;
			if (s * C + C - 1 > N)
				subvector_size = N - s * C;
			result = add_mod(mul_mod(result, powers_of_w.back(), mod), 
				reduce_subvector(&V[s * C], subvector_size, mod, powers_of_w.data()), mod);
			if (s < T)
				break;
		}
		partial_results[t].value = result;
	};
	std::vector<std::thread> threads;
	for (unsigned t = 1; t < T; ++t)
		threads.emplace_back(thread_proc, t);
	thread_proc(0);
	for (auto& thread:threads)
		thread.join();
	for (std::size_t t = 1; t < T; ++t)
		partial_results[0].value = add_mod(partial_results[0].value, mul_mod(partial_results[t].value, powers_of_w[C + t - 1], mod), mod);
	return partial_results[0].value;
}