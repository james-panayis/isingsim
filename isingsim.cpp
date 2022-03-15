/*
  Simulates the Ising model
*/

#define FMT_HEADER_ONLY 1
#include "fmt/format.h"
#include "fmt/compile.h"

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
#include <ranges>
#include <cassert>

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


// Go to the next element of a multidimensional array
template<class Space>
constexpr bool increment(std::array<std::int64_t, Space::Dimension>& indexes)
{
  for (std::size_t i = Space::Dimension - 1; i < Space::Dimension; --i)
  {
    ++indexes[i];

    if (indexes[i] >= Space::Sizes[i])
      indexes[i] = 0;
    else
      return true;
  }
  return false;
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


// Struct used to store the state; just a looping mdarray with helpful functions
template<class T, std::int64_t S0 = 0, std::int64_t... S>
struct Periodic_array
{
public:

  static_assert(S0 > 0, "Each dimension of a Periodic_array must have positive size.");

  static constexpr std::size_t Dimension{sizeof...(S) + 1};
  static constexpr std::array<std::int64_t, Dimension> Sizes{{S0, S...}};

  using Base = std::array< std::conditional_t<Dimension == 1, T, Periodic_array<T, S...>>, S0 >;
  using Sub  = Base::value_type;

  static_assert((Dimension==1) == (std::is_same_v<Sub,T>), "In a Periodic_array, the subtype should be the same as the underlying type if and only if there is only 1 dimension.");


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
  consteval std::size_t size() const
  { return base.size(); }

  // Total (recursive) number of elements
  consteval std::size_t count() const
  {
    if constexpr (Dimension == 1)
    {
      return (*this).size();
    }
    else
    {
      return base[0].count() * (*this).size();
    }
  }


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


  // Return access to view the 2D slice through the last 2 dimensions
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


  // Calculate the total megnetization in the system
  constexpr std::int64_t get_base_magnetization() const
  {
    std::int64_t magnetization{};
    std::array<std::int64_t, Dimension> indexes{};
    do
    {
      magnetization += (*this)[indexes];
    } while (increment<std::remove_reference_t<decltype(*this)>>(indexes));

    return magnetization;
  }


  // Calculate the base energy in the system
  constexpr std::int64_t get_base_energy() const
  {
    std::int64_t energy{};
    std::array<std::int64_t, Dimension> indexes{};
    do
    {
      energy += (*this).count_identical_interactions(indexes);
    } while (increment<std::remove_reference_t<decltype(*this)>>(indexes));

    assert(energy % 2 == 0);

    return energy / 2;
  }


private:

  // Store
  Base base;
};


static constexpr bool SLICE  = true;
static constexpr bool SQUISH = false;

// Create a ppm file of a state, either squishing down extra dimensions into 2 or taking a 2D slice
template<bool Slice_Instead_Of_Squish = SLICE, class Space>
void make_ppm(Space& data, const std::string filename)
{
  if constexpr (Space::Dimension <= 1)
  {
    return;
  }

  std::array<std::int64_t, Space::Dimension> indexes{};

  std::span<std::int64_t, 2> output_indexes = std::span(indexes).template last<2>();

  static constexpr std::array<std::int64_t, 2> output_sizes{Space::Sizes.end()[-2], Space::Sizes.end()[-1]};

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
  std::string tempfile = "/tmp/" + filename;
  static auto& characters = "2345689abcdefghjkmnpqrstuvwxyzABCDEFGHJKMNPQRSTUVWXYZ";
  static std::mt19937 random_generator{std::random_device{}()};
  static std::uniform_int_distribution<std::size_t> index(0, sizeof(characters) - 2);

  for (int i = 0; i < 20; i++)
    tempfile += characters[index(random_generator)];

  std::ofstream graphfile(tempfile);

  // Write data to temporary file
  graphfile << "P3 " << output[0].size() << " " << output.size() << " 255" << std::endl;
  
  static constexpr double layers = Slice_Instead_Of_Squish ? 1 : std::accumulate(begin(Space::Sizes), end(Space::Sizes)-2, 1.0, std::multiplies<double>());

  for (auto& row : output)
  {
    for (auto& bit : row)
    {
      graphfile << sc<int>(bit * 255.0 / layers) << " " << sc<int>(bit * 255.0 / layers) << " " << sc<int>(bit * 255.0 / layers) << " ";
    }
    graphfile << '\n';
  }

  graphfile.close();

  // Copy temporary file to final location
  std::filesystem::remove(filename);
  std::filesystem::copy(tempfile, filename);
  std::filesystem::remove(tempfile);
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

  // Specify dimension and size(s)
  using Space = Periodic_array<bool, 32, 128, 256>;//8, 128, 256>; //dimensions

  // Stringify size info
  std::string dimensions_string = "-" + std::to_string(Space::Dimension) + "D-";
  for (auto s : Space::Sizes)
    dimensions_string += std::to_string(s) + "x";
  dimensions_string.pop_back();

  // Set temperature
  double temperature_factor{};

  // Create arrays
  Space spaceorig;
  Space space1;
  Space space2 [[maybe_unused]];

  // Create random number machinery
  std::random_device source;
  const auto random_data = std::views::iota(std::size_t(), (std::mt19937::state_size * sizeof(typename std::mt19937::result_type) - 1) / sizeof(source()) + 1)
                           | std::views::transform([&](auto){ return source(); });
  std::seed_seq seeds(std::begin(random_data), std::end(random_data));
  std::mt19937 gen(seeds);
  std::uniform_real_distribution<double> uniform0to1(0.0, std::nextafter(1.0, std::numeric_limits<double>::max()));
  std::bernoulli_distribution fraction_per_pass(0.1); //fraction of states to examine on each pass

  // Map from number of identical adjacent spins to chance of flipping
  std::array<double, (2*Space::Dimension)+1> flip_probabilities;

  const auto initialize_probabilities = [&]
  {
    for (std::size_t i = 0; i < flip_probabilities.size(); ++i) 
    {
      flip_probabilities[i] = std::exp(-sc<double>((sc<std::int64_t>(i)-sc<std::int64_t>(Space::Dimension)) * 4)/(temperature_factor*2.0*Space::Dimension));
    }
  };


  // Perform one iteration of the multi-flip metropolis algorithm with input start_space, outputting to end_space
  const auto mfmstep [[maybe_unused]] = [&](const auto& start_space, auto& end_space) constexpr
  {
    std::array<std::int64_t, Space::Dimension> indexes{};

    // Loop through all points
    do
    {
      // Keep state
      if (!(fraction_per_pass(gen)))
      {
        end_space[indexes] = start_space[indexes];
        continue;
      }

      // Examine and possibly flip state
      end_space[indexes] = start_space[indexes] ^ ( flip_probabilities[start_space.count_identical_interactions(indexes)] > uniform0to1(gen));

    } while (increment<Space>(indexes));
  };

  // Create distributions to generate random indexes into the space
  auto index_generators = []
  {
    std::array<std::uniform_int_distribution<std::int64_t>, Space::Dimension> generators;

    for (std::size_t i = 0; i < generators.size(); ++i)
    {
      generators[i] = std::uniform_int_distribution<std::int64_t>(0, Space::Sizes[i] - 1);
    }
    return generators;
  }(); 

  // Perform one iteration of the original metripolis algorithm on space
  const auto mstep [[maybe_unused]] = [&](auto& space) constexpr
  {
    std::array<std::int64_t, Space::Dimension> indexes{};

    for (std::size_t i = 0; i < indexes.size(); ++i)
    {
      indexes[i] = index_generators[i](gen);
    }

    space[indexes] ^= (flip_probabilities[space.count_identical_interactions(indexes)] > uniform0to1(gen));
  };

  std::vector<std::vector<std::vector<std::int64_t>>> output_stats;

  static constexpr int repeats = 1;

  // Multiples of the critical value in the ferromagnetic model
  //static constexpr std::array<double, 10> temperatures{{0.1, 0.25, 0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5}};
  static constexpr std::array<double, 1> temperatures{{0.5}};

  for (auto& temperature : temperatures)
  {
    fmt::print(FMT_COMPILE("starting temperature {}\n"), temperature);
    temperature_factor = temperature;
    std::string temperature_string = std::to_string(temperature);

    initialize_probabilities();


    // Experiment repeat loop
    for (int repeat = 0; repeat < repeats; ++repeat)
    {
      fmt::print(FMT_COMPILE("starting repeat {}\n"), repeat);
    
      output_stats.resize(output_stats.size() + 1);

      // Randomize initial state
      spaceorig.randomize(std::bernoulli_distribution{0.5}, gen); //chance of up state
      space1 = spaceorig;

      make_ppm<SQUISH>(space1, "000000" + dimensions_string + ".ppm");

      // Main loop
      for (int i = 1; i <= 300 * 5; ++i)
      {
        mfmstep(space1, space2);
        mfmstep(space2, space1);

        //for (std::size_t j = 0; j < space1.count(); j++)
        //  mstep(space1);

        if (i%1==0)
        {
          output_stats.back().push_back({space1.get_base_energy(), space1.get_base_magnetization()});
          make_ppm<SQUISH>(space1, fmt::format(FMT_COMPILE("{:0>6}{}.ppm"), i, dimensions_string));
        }
      }
    }
  }


  // Print result stats to a csv

  std::string sf = "stats" + dimensions_string + ".csv";
  std::ofstream statsfile(sf);

  for (std::size_t i = 0; i < output_stats[0].size(); ++i)
  {
    for (std::size_t j = 0; j < output_stats.size(); ++j)
    {
      assert(output_stats[0].size() == output_stats[j].size());
      statsfile << output_stats[j][i][0] << ", " << output_stats[j][i][1] << ", ";
    }
    statsfile << "\n";
  }

  statsfile.close();

  return EXIT_SUCCESS;
}



















