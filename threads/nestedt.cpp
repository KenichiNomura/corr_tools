#include <iostream>
#include <thread>
#include <chrono>
#include <vector>
#include <mutex>
#include <omp.h>

std::mutex mtx; 

void sleep_for(int num_threads)
{
	int num_steps = 10;
	double step = 1.0 / (double)num_steps;


	std::vector<double> sum;
	sum.resize(num_threads,0.0);

	auto start = std::chrono::high_resolution_clock::now();

	#pragma omp parallel num_threads(num_threads)
	{
		int tid = omp_get_thread_num();

		// do some work, like GPU offloading
		for(int i = tid; i < num_steps; i += num_threads)
			std::this_thread::sleep_for(std::chrono::seconds(1));

	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end-start;

	for(int tid = 1; tid < num_threads; tid++) sum[0] += sum[tid];

	std::lock_guard<std::mutex> guard(mtx);
	std::cout << "pi = " << step * sum[0] << " by " << num_threads << " threads " << elapsed.count() << std::endl;;

}


struct NestedThreads
{

	int num_threads;
	std::vector<std::thread> threads;

	NestedThreads(int _num_threads=1) : num_threads(_num_threads) {}

	void run()
	{
		for(int ithread=0; ithread<num_threads; ithread++)
			threads.push_back( std::thread(sleep_for, ithread+1) );

		for (auto &th : threads) 
			th.join();
	}

};

int main(int argc, char *argv[])
{
	NestedThreads l(std::atoi(argv[1]));

	l.run();
}
