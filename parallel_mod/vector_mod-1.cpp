#include "vector_mod.h"
#include "mod_ops.h"
#include "num_threads.h"
#include <vector>
#include <thread>
#include <cassert>

#if defined(__GNUC__) && __GNUC__ < 10
namespace std
{
	constexpr std::size_t hardware_constructive_interference_size = 64u;
	constexpr std::size_t hardware_destructive_interference_size = 64u; //See False sharing
}
#endif

constexpr std::size_t C = std::hardware_constructive_interference_size;
const std::size_t T = get_num_threads();

static std::vector<IntegerWord> get_powers_of_w(IntegerWord mod)
{
	//IntegerWord A = 1;
	std::vector<IntegerWord> table;
	table.emplace_back(1);
	for (std::size_t i = 1; i <= C; ++i)
		table.emplace_back(times_word(table.back(), mod));
	for (std::size_t t = 2; t <= T; ++t)
		table.emplace_back(mul_mod(table.back(), table[C]. mod));
	return table;
}

static IntegerWord reduce_subvector(const IntegerWord* subvector, std::size_t count, IntegerWord mod)
{
	verify(count <= C);
	IntegerWord A = 0;
	while (count--)
	{
		A = times_word(A, mod);
		A = add_mod(A, subvector[count], mod);
	}
	//See Horner method
	return A;
}

IntegerWord vector_mod(const IntegerWord* V, std::size_t N, IntegerWord mod)
{
	const auto powers_of_w = get_powers_of_w(mod);
	std::vector<std::thread> threads;
	struct element_t
	{
		alignas(std::hardware_destructive_interference_size) IntegerWord value;
	};
	std::vector<element_t> partial_results(T);
	auto thread_proc = [&powers_of_w, &partial_results](std::size_t t)
	{
	};
	for (std::size_t t = 1; t < T; ++t)
		threads.emplace_back(thread_proc, t);
	thread_proc(0);
	for (auto& thread:threads)
		thread.join();
	for (std::size_t i = 1; i < T; ++i)
		partial_results[0].value = add_mod(partial_results[0].value, mul_mod(partial_results[i],  powers_of_w[C + i], mod), mod);
	return partial_results[0].value;
}
