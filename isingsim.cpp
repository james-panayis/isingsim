/*
  Simulates the Ising model
*/

#define FMT_HEADER_ONLY 1
#include "fmt/format.h"

#include <iostream>
#include <fstream>
#include <string>
#include <charconv>
#include <array>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <cmath>
#include <vector>
#include <filesystem>
#include <random>
#include <span>

// Shorter static_cast
template<class T>
constexpr auto sc(auto a)
{
  return static_cast<T>(a);
}

//get the Nth element of a parameter pack
template<std::size_t I, class... Ts>
constexpr auto get(Ts&&... ts)
{
  return std::get<I>(std::forward_as_tuple(ts...));
}


// Signed modulo
constexpr std::int64_t mod(const std::int64_t value, const std::int64_t max)
{
  return (value % max + max) % max;
}

// Clamp
constexpr std::int64_t clm(const std::int64_t value, const std::int64_t max)
{
  return std::clamp(value, std::int64_t{0}, max);
}

// Identity
constexpr std::int64_t idt(const std::int64_t value, const std::int64_t max [[maybe_unused]])
{
  return value;
}

// Create Periodic_array struct

// Forward declaration of internal Periodic_array implementaion
template<class T, std::size_t Dimension, std::int64_t... Sizes>
struct Periodic_array_internal;

// Return an object of the type that a Periodic_array should contain
// This function should never be called
template<class T, std::size_t Dimension, std::int64_t Size0, std::int64_t... Other_sizes>
consteval auto make_periodic_array_internal_base()
{
  //std::unreachable(); //C++23, not yet implemented

  static_assert(Dimension == sizeof...(Other_sizes) + 1, "Dimensional mismatch in make_periodic_array_internal_base_type");

  if constexpr (Dimension != 1)
  {
    return std::array<Periodic_array_internal<T, Dimension - 1, Other_sizes...>, Size0>();
  }
  else
  {
    return std::array<T, Size0>();
  }
}

// Use return type of above function
template<class T, std::size_t Dimension, std::int64_t... Sizes>
using Periodic_array_internal_base = decltype(make_periodic_array_internal_base<T, Dimension, Sizes...>());

// Internal implementation of Periodic_array
template<class T, std::size_t N, std::int64_t... S>
struct Periodic_array_internal
{
  static constexpr std::size_t Dimension = N;
  static constexpr std::array<std::int64_t, Dimension> Sizes{{S...}};

  using Base = Periodic_array_internal_base<T, Dimension, S...>;
  using Sub  = Base::value_type;

  // Subtype should be underlying type if and only if we are 1-dimensional
  static_assert((Dimension==1) == (std::is_same_v<Sub, T>));


  // Subscript operators
  constexpr Sub& operator[](std::int64_t index)
  { return base[mod(index, base.size())]; }

  constexpr const Sub& operator[](std::int64_t index) const
  { return base[mod(index, base.size())]; }

  template<std::integral... Indexes>
  constexpr auto& operator[](std::int64_t index, Indexes... other_indexes)
  { return (*this)[index][other_indexes...]; }

  template<std::integral... Indexes>
  constexpr const auto& operator[](std::int64_t index, Indexes... other_indexes) const
  { return (*this)[index][other_indexes...]; }

  constexpr auto& operator[](std::span<std::int64_t, Dimension> indexes)
  { 
    if constexpr (Dimension == 1)
    {
      return (*this)[indexes[0]];
    }
    else
    {
      return (*this)[indexes[0]][indexes.template last<Dimension-1>()];
    }
  }

  constexpr const auto& operator[](std::span<std::int64_t, Dimension> indexes) const
  { 
    if constexpr (Dimension == 1)
    {
      return (*this)[indexes[0]];
    }
    else
    {
      return (*this)[indexes[0]][indexes.template last<Dimension-1>()];
    }
  }


  // Iterators
  constexpr auto begin()
  { return base.begin(); }

  constexpr auto end()
  { return base.end(); }

  constexpr const auto begin() const
  { return base.begin(); }

  constexpr const auto end() const
  { return base.end(); }


  // Size
  constexpr std::size_t size() const
  { return base.size(); }


  // Fill with random data
  constexpr void randomize(auto distribution, std::mt19937& generator)
  {
    if constexpr (Dimension == 1)
    {
      for (T& value : base)
        value = distribution(generator);
    }
    else
    {
      for (Sub& subspace : base)
        subspace.randomize(distribution, generator);
    }
  }


  // Calculate how many neighbours (interacting points) of a point are in the same state/spin as the point
  constexpr std::int64_t count_identical_interactions(std::span<std::int64_t, Dimension> indexes) const
  {
    if constexpr (Dimension == 1)
    {
      return  ((*this)[indexes[0]-1] == (*this)[indexes[0]])
             +((*this)[indexes[0]+1] == (*this)[indexes[0]]);
    }
    else
    {
      auto subindexes = indexes.template last<Dimension - 1>();

      return  ((*this)[indexes[0]-1][subindexes] == (*this)[indexes[0]][subindexes])
             +((*this)[indexes[0]+1][subindexes] == (*this)[indexes[0]][subindexes])
             +(*this)[indexes[0]].count_identical_interactions(subindexes);
    }
  }


  // Return the 2D slice through the last 2 dimensions
  constexpr const auto& slice2D() const
  {
    if constexpr (Dimension == 2)
    {
      return (*this);
    }
    else
    {
      return base[0].slice2D();
    }
  }


  // Store
  Base base;
};

// Exposed Periodic_array (does not need Dimension template parameter, unlike Periodic_array_internal)
template<class T, std::int64_t... Sizes>
using Periodic_array = Periodic_array_internal<T, sizeof...(Sizes), Sizes...>;


// Go to the next element of a multidimensional array
template<class Space>
constexpr bool increment(std::array<std::int64_t, Space::Dimension>& indexes)
{
  for (std::size_t i = 0; i < Space::Dimension; ++i)
  {
    ++indexes[i];

    if (indexes[i] >= Space::Sizes[i])
      indexes[i] = 0;
    else
      return true;
  }
  return false;
}


constexpr bool SLICE  = true;
constexpr bool SQUISH = false;

// Create a ppm file of a state, either squishing down extra dimensions into 2 or taking a 2D slice
template<bool Slice_Instead_Of_Squish = SLICE, class Space>
void make_ppm(Space& data, const std::string filename)
{
  std::array<std::int64_t, Space::Dimension> indexes{};

  std::span<std::int64_t, 2> output_indexes = std::span(indexes).template last<2>();

  constexpr std::array<std::int64_t, 2> output_sizes{Space::Sizes.end()[-2], Space::Sizes.end()[-1]};

  std::array<std::array<double, output_sizes[1]>, output_sizes[0]> output{};

  // Loop through all points
  do
  {
    if constexpr (Slice_Instead_Of_Squish)
    {
      output[output_indexes[0]][output_indexes[1]] = data.slice2D()[output_indexes];
    }
    else
    {
      output[output_indexes[0]][output_indexes[1]] += data[indexes];
    }
  } while (increment<Space>(indexes));

  // Create temporary file
  std::string tmp = "graphtest";
  //std::stringstream tomakeuid(loc);
  //while(std::getline(tomakeuid, tmp, '/')) {}

  std::ofstream graphfile("/tmp/" + tmp);

  std::cout << "P3 " << output[0].size() << " " << output.size() << " 255" << std::endl;
  graphfile << "P3 " << output[0].size() << " " << output.size() << " 255" << std::endl;
  
  constexpr double layers = Slice_Instead_Of_Squish ? 1 : std::accumulate(begin(Space::Sizes), end(Space::Sizes)-2, 1.0, std::multiplies<double>());

  for (auto& row : output)
  {
    for (auto& bit : row)
    {
      graphfile << sc<int>(bit * 255.0 / layers) << " " << sc<int>(bit * 255.0 / layers) << " " << sc<int>(bit * 255.0 / layers) << " ";
    }
    graphfile << '\n';
  }

  std::cout << "done generating ppm" << std::endl;

  graphfile.close();

  // Copy temporary file to final location
  std::filesystem::remove(filename);
  std::filesystem::copy("/tmp/" + tmp, filename);
  std::filesystem::remove("/tmp/" + tmp);
}


int main(int argc, char* argv[])
{
  // Check number of arguments passed
	if (argc != 1)
	{
		std::cout<<argc << "usage: " << argv[0] << '\n';
		return EXIT_FAILURE;
	}

	/*int return_value = EXIT_SUCCESS;

	if (argc > 2)
	{
		for (int i = 1; i < argc; i++)
		{
			char* argv2[] = {argv[0], argv[i]};

			std::cout << argv2[0] << " " << argv2[1] << std::endl;
			
			#pragma GCC diagnostic push
   		#pragma GCC diagnostic ignored "-Wpedantic" //disable warning for (non-iso-compliant) calling of main

			return_value = (main(2, argv2) == return_value) ? return_value : EXIT_FAILURE;

			#pragma GCC diagnostic pop
		}
		return return_value;
	}*/

  // Specify dimension and sizes
  using Space = Periodic_array<bool, 16, 128, 128>;//8, 128, 256>; //dimensions

  // Create arrays
  Space space1;
  Space space2 [[maybe_unused]];

  // Create random number machinery
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> u(0.0, std::nextafter(1.0, std::numeric_limits<double>::max()));
  std::bernoulli_distribution rb(0.12); //fraction of states to examine on each pass

  // Randomize initial state
  space1.randomize(std::bernoulli_distribution{0.5}, gen); //chance of up state

  // Map from number of identical adjacent spins to chance of flipping
  constexpr auto flip_probabilities = []
  {
    std::array<double, (2*Space::Dimension)+1> probability_table{};

    for (std::size_t i = 0; i < probability_table.size(); i++) 
    {
      probability_table[i] = std::exp(-sc<double>((sc<std::int64_t>(i)-sc<std::int64_t>(Space::Dimension)) * 4)/1.0);
    }
    return probability_table;
  }();

  // Perform one iteration with input start_space, outputting to end_space
  const auto step = [&](auto& start_space, auto& end_space) constexpr
  {
    /////// NEW CODE ///////
    
    std::array<std::int64_t, Space::Dimension> indexes{};

    // Loop through all points
    do
    {
      // Keep state
      if (!(rb(gen)))
      {
        end_space[indexes] = start_space[indexes];
        continue;
      }

      // Examine and possibly flip state
      end_space[indexes] = start_space[indexes] ^ ( flip_probabilities[start_space.count_identical_interactions(indexes)] > u(gen));

    } while (increment<Space>(indexes));

    /////// OLD CODE ///////

    /*for (std::size_t d = 0; d < start_space.size(); d++)
    {
      for (std::size_t r = 0; r < start_space[0].size(); r++)
      {
        for (std::size_t c = 0; c < start_space[0][0].size(); c++)
        {
          if (!(rb(gen)))
          {
            end_space[d,r,c] = start_space[d,r,c];
            continue;
          }

          const double Ediff = -(
                                  (start_space[d,r,c] == start_space[d  ,r,  c-1])
                                + (start_space[d,r,c] == start_space[d  ,r,  c+1])
                                + (start_space[d,r,c] == start_space[d  ,r-1,c  ])
                                + (start_space[d,r,c] == start_space[d  ,r+1,c  ])
                                + (start_space[d,r,c] == start_space[d+1,r  ,c  ])
                                + (start_space[d,r,c] == start_space[d-1,r  ,c  ])
                                //+ (start_space[r,c] == start_space[r+1,c+1])
                                //+ (start_space[r,c] == start_space[r-1,c-1])
                                //+ (start_space[r,c] == start_space[r+1,c-1])
                                //+ (start_space[r,c] == start_space[r-1,c+1])
                                //+ (start_space[r,c] == start_space[r-1,c-1])
                                //- 2) * -4;
                                - 3) * -4;


          end_space[d,r,c] = start_space[d,r,c] 
                           ^ ( ((Ediff<0) || ((memo_exp(-Ediff/1.0) > u(gen)))) 
                               //&& rb(gen)
                                              );
                               //&& (0.8225 > u(gen)) );

        }
      }
    }*/
  };

  make_ppm<SQUISH>(space1, "000000.ppm");

  // Main loop
  for (int i = 1; i <= 900; i++)
  {
    step(space1, space2);

    //make_ppm(space2, fmt::format("{:0>6}", i)+"a.ppm");

    step(space2, space1);

    if (i%3==0) make_ppm<SQUISH>(space1, fmt::format("{:0>6}", i)+".ppm");
  }

  return EXIT_SUCCESS;
}




















