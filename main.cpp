#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <iterator>



typedef std::tuple<double, double, double> coord;

std::string exec (const std::string& str);

std::vector<coord> database_read (const std::string& name);

std::vector<std::pair<double, double>> analyzer (const std::string& name,
                                                 std::vector<std::pair<double, double>>& distribution_function);

void plots (const std::string& name, int& time_step);


// The code bellow using for transformation str to num.
template<typename T>
T fromString (const std::string& s) {
    std::istringstream iss(s);
    T res;
    iss >> res;
    return res;
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


template <typename T>
void data_file_creation (const std::string& name, const T& xx) {
    std::ofstream fout;
    fout.open(name);
    for (int i = 0; i < xx.size(); ++i)
        if (std::isfinite(xx.at(i).first) && std::isfinite(xx.at(i).second))
            fout << xx.at(i).first << '\t' << xx.at(i).second << std::endl;
    fout.close();
}


int main () {
    std::string distributions_files_path = std::move(exec("rm -rf distributions && mkdir distributions && cd distributions && echo $PWD"));
    std::string analyzed_data_files_name = distributions_files_path + "/Distribution";
    std::string source_files_name = "/home/alexander/CLionProjects/StatAnalyzer/cmake-build-debug/velocities/Ar_velocities";
    int frames_count = fromString<int>(exec("cd /home/alexander/CLionProjects/StatAnalyzer/cmake-build-debug/velocities/ && ls | wc -l"));
//#pragma omp parallel
    {

        for (int i = 0; i < frames_count; ++i) {

            std::vector<std::pair<double, double>> distribution_function;
            std::vector<std::pair<double, double>> data = std::move(analyzer(source_files_name + '.' + toString(i),
                                                                             distribution_function));
            std::string name = analyzed_data_files_name + '.' + toString(i);
            data_file_creation(name, data);
            data_file_creation(name + '.' + "analytic", distribution_function);
            plots(analyzed_data_files_name, i);
        }
    }
}


bool is_equal (double a, double b) {
    return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}


template<typename T, size_t... Is>
double vector_length_impl (T const& t, std::index_sequence<Is...>) {
    return std::sqrt((std::pow(std::get<Is>(t), 2) + ...));
}

template <class Tuple>
double vector_length (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return vector_length_impl(t, std::make_index_sequence<size>{});
}


int group_definition (std::vector<double>& groups_borders, double& velocity) {
    for (int i = groups_borders.size(); i > 0; --i)
        if (velocity <= groups_borders[i] && velocity >= groups_borders[i-1])
            return i-1;
}



std::vector<double> SimpleSplineMoments(std::vector<double>& f, double& h){
    std::vector<double> m;
    int n = f.size()-1;
    m.push_back((3*f[0] - 4*f[1] + f[2])/2.0/h);
    for(int i = 1; i < n; ++i)
        m.push_back((f[i+1] - f[i-1])/2.0/h);
    m.push_back((-3*f[n] + 4*f[n-1] - f[n-2])/2.0/h);
    return m;
}

std::vector<double> CubicSpline(std::vector<double>& x, std::vector<double>& f, std::vector<double>& xx){
    std::vector<double> yy;
    int i = 0;
    double h = x[1] - x[0];
    std::vector<double> m = SimpleSplineMoments(f, h);
    for(int j = 0; j < x.size()-1; ++j)
        while(xx[i] >= x[j] && xx[i] <= x[j+1]){
            double buf1 = x[j+1] - xx[i];
            double buf2 = xx[i] - x[j];
            yy.push_back(std::pow(buf1, 2)*(2*buf2+h)*f[j]/std::pow(h, 3)+
                         std::pow(buf2, 2)*(2*buf1+h)*f[j+1]/std::pow(h, 3)+
                         std::pow(buf1, 2)*buf2*m[j]/std::pow(h, 2)-
                         std::pow(buf2, 2)*buf1*m[j+1]/std::pow(h, 2));
            ++i;
        }
    return yy;
}

const double pi = 3.14159265359;
const double m = 6.6335e-23;
const double k_B = 1.380649e-16;
const double T = 300;

double maxwell_distribution (double& v) {
    return 4.0*pi*std::pow(v, 2) * std::pow((m / (2.0 * pi * k_B * T)), 3.0/2.0)
                                    * std::exp(- m * std::pow(v, 2) / (2.0 * k_B * T));
}


template<typename T, size_t... Is>
bool is_finite_tuple_impl (T const& t, std::index_sequence<Is...>) {
    return (std::isfinite(std::get<Is>(t)) & ...);
}

template <class Tuple>
bool is_finite_tuple (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return is_finite_tuple_impl(t, std::make_index_sequence<size>{});
}



std::vector<std::pair<double, double>> analyzer (const std::string& name,
                                                 std::vector<std::pair<double, double>>& distribution_function) {

    std::vector<coord> data = database_read(name);
    int number_of_particles = data.size();

    std::vector<double> abs_velocities;


    for (int i = 0; i < number_of_particles; ++i) {
        if (is_finite_tuple(data[i]))
            abs_velocities.emplace_back(vector_length(data[i]));
    }
    std::sort(abs_velocities.begin(), abs_velocities.end());


    double V_max = *max_element(abs_velocities.begin(), abs_velocities.end());
    double V_min = *min_element(abs_velocities.begin(), abs_velocities.end());


    double spline_step = 1000;

    std::vector<double> vv (spline_step);

    double V_0 = V_min;
    double V_spline_step = (V_max - V_min) / spline_step;
    std::generate(vv.begin(), vv.end(), [&] {V_0 += V_spline_step; return V_0;});

    std::vector<double> f;
    for (double & abs_velocity : abs_velocities)
        f.emplace_back(maxwell_distribution(abs_velocity));

    std::vector<double> analytic_distribution = std::move(CubicSpline(abs_velocities, f, vv));

    for (int i = 0; i < analytic_distribution.size(); ++i)
        distribution_function.emplace_back(std::make_pair(vv[i], analytic_distribution[i]));



    if (is_equal(V_min, V_max)) {
        std::vector<std::pair<double, double>> distribution;
        distribution.emplace_back(std::make_pair(1, V_min));
    } else {


        //Cubic spline for analytic form of Maxwell Distribution


        double number_of_groups = number_of_particles / 10.0; // This param is double because will be used later as double.

        std::vector<double> groups_borders (number_of_groups+1);

        double V_step = (V_max - V_min) / number_of_groups;

        V_0 = V_min;
        std::generate(groups_borders.begin(), groups_borders.end(), [&] {V_0 += V_step; return V_0;});

        std::vector<int> numbers (number_of_particles);
        for (int i = 0; i < number_of_particles; ++i)
            numbers[i] = group_definition(groups_borders, abs_velocities[i]);

        std::vector<std::pair<double, double>> distribution (number_of_groups);

        for (int i = 0; i < number_of_particles; ++i)
            for (int j = 0; j < groups_borders.size()-1; ++j)
                if (numbers[i] == j) {
                    distribution[j].first += 1.0 / number_of_particles;
                    distribution[j].second = groups_borders[j];
                }

        std::sort(distribution.begin(), distribution.end(), [] (auto& left, auto& right)
        {return left.second < right.second;});

        return distribution;
    }
}


namespace std {
    istream& operator >> (istream& in, coord& data) {
        double first, second, third;
        in >> first >> second >> third;
        data = {first, second, third};
        return in;
    }

    ostream& operator << (ostream& out, const coord& data) {
        auto [first, second, third] = data;
        out << first << ' ' << second << ' ' << third << ' ';
        return out;
    }
}


std::vector<coord> database_read (const std::string& name) {
    std::ifstream inFile(name);
    std::vector<coord> tuples_vector;
    copy(std::istream_iterator<coord> {inFile}, std::istream_iterator<coord> {},
         back_inserter(tuples_vector));
    //copy(tuples_vector.begin(), tuples_vector.end(), std::ostream_iterator<coord>(std::cout, "\n"));
    return tuples_vector;
}


//The function returns the terminal ans. Input - string for term.
std::string exec (const std::string& str) {
    const char* cmd = str.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    result = result.substr(0, result.length()-1);
    return result;
}

void plots (const std::string& name, int& time_step) {
    FILE *gp = popen("gnuplot  -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 700, 700",
                                      "set output \'" + name + toString(time_step) + ".jpg\'",
                                      "set title \'Time step: " + toString(time_step) + "\'",
                                      "set grid xtics ytics",
                                      "set style fill solid",
                                      "set yrange [0:1]",
                                      "set key off",
                                      "set border 4095",
                                      "plot \'" + name + "." + toString(time_step) + "\' using 2:1 w boxes,\
                                      \'" + name + "." + toString(time_step) + "." + "analytic" + "\' using 1:2 with lines"};
    for (const auto& it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}
