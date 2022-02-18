/*
  Simulates the Ising model
*/

#define FMT_HEADER_ONLY 1
#include "fmt/format.h"

#include "memo/memo.hpp"

#include <cxxabi.h>

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

/*// Define mdarray

// This should never be called; only exists to generate return type
template<class T, std::size_t Dimension, std::array<std::int64_t, Dimension> Sizes>
consteval auto make_mdarray() {
  throw "calling function make_mdarray";
  if constexpr (Dimension != 1)
  {
    constexpr std::array<std::int64_t, Dimension - 1> new_indexes = [&]()
    {
      std::array<std::int64_t, Dimension - 1> n;
      for (std::size_t i = 0; i < Dimension - 1; i++)
        n[i] = Sizes[i+1];
      return n;
    }();
    return std::array<decltype(make_mdarray<T, Dimension - 1, new_indexes>()), Sizes[0]>();
  }
  else
  {
    return std::array<T, Sizes[0]>();
  }
}

template<class T, std::size_t Dimension, std::array<std::int64_t, Dimension> Sizes>
using mdarray = decltype(make_mdarray<T, Dimension, Sizes>());*/


// Signed modulo
constexpr std::int64_t mod(const std::int64_t value, const std::int64_t max)
{
  return (value % max + max) % max;
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


  // Calculate energy change caused by flipping a value
  constexpr std::int64_t flip_energy(std::span<std::int64_t, Dimension> indexes) const
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
             +(*this)[indexes[0]].flip_energy(subindexes);
    }
  }


  // Store
  Base base;
};

// Exposed Periodic_array (does not need Dimension template parameter, unlike Periodic_array_internal)
template<class T, std::int64_t... Sizes>
using Periodic_array = Periodic_array_internal<T, sizeof...(Sizes), Sizes...>;


template<std::size_t n>
using Vec = std::array<float, n>;

template<std::size_t n>
Vec<n> operator-(Vec<n> a, Vec<n> b)
{
	Vec<n> out{};

	out[0] = a[0] - b[0];
	out[1] = a[1] - b[1];
	out[2] = a[2] - b[2];

	return out;
}

template<std::size_t n>
Vec<n> operator+(Vec<n> a, Vec<n> b)
{
	Vec<n> out{};

	out[0] = a[0] + b[0];
	out[1] = a[1] + b[1];
	out[2] = a[2] + b[2];

	return out;
}

template<std::size_t n>
Vec<n> operator/(Vec<n> a, float b)
{
	Vec<n> out{};

	out[0] = a[0] / b;
	out[1] = a[1] / b;
	out[2] = a[2] / b;

	return out;
}

template<std::size_t n>
Vec<n> operator*(Vec<n> a, float b)
{
	Vec<n> out{};

	out[0] = a[0] * b;
	out[1] = a[1] * b;
	out[2] = a[2] * b;

	return out;
}

template<std::size_t n>
Vec<n> operator*(float b, Vec<n> a)
{
	return a * b;
}

//returns the square magnitude of a Vec
template<std::size_t n>
float dot(Vec<n> v, Vec<n> w)
{
	float val{0};
	
	for (std::size_t i = 0; i < n; i++)
		val += v[i] * w[i];
	
	return val;
}

//returns the square magnitude of a Vec
template<std::size_t n>
inline float getsqmag(Vec<n> v)
{
	return dot(v, v);
}

//returns the magnitude of a Vec
template<std::size_t n>
inline float getmag(Vec<n> v)
{
	return std::sqrt(getsqmag(v));
}


//returns a Vec<3> perpendicular to the input Vec<3>
Vec<3> getperp(Vec<3> v)
{
	Vec<3> out;

	out[0] = 1.0f;
	out[1] = 1.0f;
	out[2] = - (v[0] + v[1]) / (v[2] + 0.00001f);

	out = out / getmag(out);

	return out;
}


// Return minimum distance between line segment vw and point p
float dist(Vec<2> v, Vec<2> w, Vec<2> p) {
  const float l2 = getsqmag(v - w);  // i.e. |w-v|^2 -  avoid a sqrt
  if (l2 == 0.0f) return getmag(p - v);   // v == w case
  // Consider the line extending the segment, parameterized as v + t (w - v)
  // Project point p onto the line, which is where t = [(p-v) . (w-v)] / |w-v|^2
  // Clamp t from [0,1] to handle points outside the segment vw
  const float t = std::clamp(dot(p - v, w - v) / l2, 0.0f, 1.0f);
  const Vec<2> projection = v + t * (w - v);  // Projection falls on the segment
  return getmag(p - projection);
}


using Pixel = std::array<std::uint32_t, 3>;

Pixel combine(Pixel p, Pixel q, float f)
{
	return Pixel{
		static_cast<std::uint32_t>((1-f) * static_cast<float>(p[0]) + f * static_cast<float>(q[0])),
		static_cast<std::uint32_t>((1-f) * static_cast<float>(p[1]) + f * static_cast<float>(q[1])),
		static_cast<std::uint32_t>((1-f) * static_cast<float>(p[2]) + f * static_cast<float>(q[2]))
	};
}

using Image = std::vector<std::vector<Pixel>>;

template<class F>
concept floatable = std::convertible_to<F, float>;

//template<std::convertible_to<float>... F>// requires (std::convertible_to<F, float>)
void drawlinei(const Image& data, Vec<2> p1, Vec<2> p2, floatable auto thickness, Pixel color, floatable auto opacity)
{
	for (std::size_t i = 0; i < data.size(); i++)
		for (std::size_t j = 0; j < data[i].size(); j++)
		{
			//float dist = std::numeric_limits<float>::max();

			Vec<2> pc = {
				(static_cast<float>(j) + 0.5f) / static_cast<float>(data[i].size()),
				(static_cast<float>(i) + 0.5f) / static_cast<float>(data.size())
			};
			
			/*for (float d = 0; d <= 1; d += 1.0f/16.0f)
			{
				float x = (1.0f-d)*p1[0] + d*p2[0];
				float y = (1.0f-d)*p1[1] + d*p2[1];
				
				dist = std::min(dist, (x-pc[0])*(x-pc[0]) + (y-pc[1])*(y-pc[1]));
			}*/

			float d = dist(p1, p2, pc);


			if (d <= thickness / 2)
				data[i][j] = combine(data[i][j], color, opacity);
			else
				if (d <= thickness)
					data[i][j] = combine(data[i][j], color, opacity * 2.0f * (thickness - d) / thickness);

			//data[i][j] = combine(data[i][j], Pixel{0, 0, static_cast<uint32_t>(dist * 255.0f)}, 0.5);
		}


	return; 
}

inline void drawline(const Image& data, Vec<2> p1, Vec<2> p2, floatable auto thickness, Pixel color, floatable auto opacity)
{
	drawlinei(data, p1, p2, static_cast<float>(thickness), color, static_cast<float>(opacity));
}






// returns the average of the locations of all the 'ATOM's in a pdb file
Vec<3> get_COM(const std::string filename)
{
	std::ifstream infile(filename);

	Vec<3> average{};
	int count = 0;

	for (std::string line; std::getline(infile, line); )
	{
		if (!line.starts_with("ATOM"))
			continue;

		Vec<3> pos;

		std::string xs = line.substr(30, 8);
		std::string ys = line.substr(38, 8);
		std::string zs = line.substr(46, 8);

		auto errx = sscanf(xs.c_str(), "%f", &pos[0]);
		auto erry = sscanf(ys.c_str(), "%f", &pos[1]);
		auto errz = sscanf(zs.c_str(), "%f", &pos[2]);

		if ((errx == EOF) || (erry == EOF) || (errz == EOF))
		{
			std::cout << "ERROR: coordinate missing in line: " << line << std::endl;
			continue;
		}

		average = average + pos;

		count++;
	}

	infile.close();

	return average / static_cast<float>(count);
}


void make_ppm(auto& data, const std::string filename)
{
  std::string tmp = "graphtest";
  //std::stringstream tomakeuid(loc);
  //while(std::getline(tomakeuid, tmp, '/')) {}

  std::ofstream graphfile("/tmp/" + tmp);

  std::uint32_t stretch = 1;

  //std::cout << "coords: " << data.size() << "  " << data[0].size() << "  " << data[1].size() << std::endl;
  std::cout << "P3 " << data[0].size() * stretch << " " << data.size() * stretch << " 255" << std::endl;
  graphfile << "P3 " << data[0].size() * stretch << " " << data.size() * stretch << " 255" << std::endl;
  //std::cout << "P3 " << 1000 * stretch << " " << 1000 * stretch << " 255" << std::endl;
  //graphfile << "P3 " << 1000 * stretch << " " << 1000 * stretch << " 255" << std::endl;

  for (auto& row : data)
  {
    for (std::uint32_t i = 0; i < stretch; i++)
    {
      for (auto& bit : row)
      {
        for (std::uint32_t j = 0; j < stretch; j++)
        {
          graphfile << sc<int>(bit) * 255 << " " << sc<int>(bit) * 255 << " " << sc<int>(bit) * 255 << " ";
        }
        //if (i==0) std::cout << bit << " ";
        //if (i==0) fmt::print("{} ", sc<int>(bit));
      }
      graphfile << '\n';
    }
    //std::cout << '\n';
  }

  std::cout << "done" << std::endl;

  graphfile.close();

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
  using Space = Periodic_array<bool, 8, 128, 256>; //dimensions

  // Create arrays
  Space space1;
  Space space2 [[maybe_unused]];

  // Create random number generators
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> u(0.0, std::nextafter(1.0, std::numeric_limits<double>::max()));
  std::bernoulli_distribution rb(0.1); //fraction of states to examine on each pass

  // Randomize initial state
  space1.randomize(std::bernoulli_distribution{0.5}, gen); //chance of up state

  // Print initial frame
  make_ppm(space1.base[0], "000000.ppm");

  // Memoize exponential function (minor performance hack)
  auto memo_exp = memo::memoize( [](const double x) constexpr { return std::exp(x); } );

  const auto step = [&](auto& start_space, auto& end_space) constexpr
  {
    /////// NEW CODE ///////
    //
    std::array<std::int64_t, Space::Dimension> indexes = {};

    const auto increment = [&]() constexpr -> bool
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
    };

    do
    {
      if (!(rb(gen)))
      {
        end_space[indexes] = start_space[indexes];
        continue;
      }

      const double Ediff = sc<double>(-(start_space.flip_energy(indexes) - 3) * -4);

      end_space[indexes] = start_space[indexes] ^ ( ((Ediff<0) || ((memo_exp(-Ediff/1.0) > u(gen)))));

    } while (increment());

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


  // Main loop
  for (int i = 1; i < 3000; i++)
  {
    step(space1, space2);

    //make_ppm(space2.base[0], fmt::format("{:0>6}", i)+"a.ppm");

    step(space2, space1);

    if (i%5==0) make_ppm(space1.base[0], fmt::format("{:0>6}", i)+".ppm");
  }

  return EXIT_SUCCESS;
}




















