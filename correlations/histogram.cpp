#include <boost/histogram.hpp>
#include <boost/format.hpp> // used here for printing
#include <iostream>

int main() {

    using namespace boost::histogram;

    // make 1d histogram with 4 regular bins from 0 to 2
    auto h = make_histogram( axis::regular<>(4, 0.0, 2.0) );

	std::cout << typeid(h).name() << std::endl;
 
    // push some values into the histogram
    for (auto&& value : { 0.4, 1.1, 0.3, 1.7, 10. })
      h(value);

	std::cout << "size " << h.axis().size() << std::endl;
    for(int i=0; i<h.axis().size(); i++)
	std::cout << h[i] << std::endl;

    // iterate over bins
    for (auto&& x : indexed(h)) {
      std::cout << boost::format("bin %i [ %.1f, %.1f ): %i\n")
        % x.index() % x.bin().lower() % x.bin().upper() % *x;
    }

    std::cout << std::flush;

    /* program output:

    bin 0 [ 0.0, 0.5 ): 2
    bin 1 [ 0.5, 1.0 ): 0
    bin 2 [ 1.0, 1.5 ): 1
    bin 3 [ 1.5, 2.0 ): 1
    */
}
