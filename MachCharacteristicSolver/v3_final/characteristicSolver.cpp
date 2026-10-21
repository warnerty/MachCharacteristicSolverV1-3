// g++ -O3 -std=c++17 -fopenmp characteristicSolver.cpp -o characteristicSolver.exe
// .\characteristicSolver.exe
// Get-Item .\characteristics.csv

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <omp.h>

// Developed by **Tyler Warner**

constexpr double PI = 3.14159265358979323846;

struct FanState
{
    double Mach;
    double theta;
    double nu;
    double mu;
    double alpha;
    double pressure;
    double temperature;
};

// convert degrees to radians
double radians(double angle)
{return angle * PI / 180.0;}

// convert radians to degrees
double degrees(double angle)
{return angle * 180.0 / PI;}

// calculate the prandtl-meyer angle
double prandtlMeyer(double Mach, double gamma)
{
    if (!std::isfinite(Mach) || Mach < 1.0 ||
        !std::isfinite(gamma) || gamma <= 1.0)
    {throw std::invalid_argument("Mach must be >= 1 and gamma must be > 1.");}

    const double a = (gamma + 1.0) / (gamma - 1.0);
    const double b = Mach * Mach - 1.0;

    return std::sqrt(a) * std::atan(std::sqrt(b / a))
         - std::atan(std::sqrt(b));
}

// find mach using bisection
double machFromNu(double nu, double gamma)
{
    if (!std::isfinite(gamma) || gamma <= 1.0)
    {throw std::invalid_argument("Gamma must be > 1.");}

    const double nuMax =
        (PI / 2.0) *
        (std::sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0);

    if (!std::isfinite(nu) || nu < 0.0 || nu >= nuMax)
    {throw std::invalid_argument("Prandtl-Meyer angle is outside the valid range.");}

    if (nu == 0.0)
    {return 1.0;}

    double lower = 1.0;
    double upper = 2.0;

    // find an upper bound
    while (prandtlMeyer(upper, gamma) < nu)
    {
        upper *= 2.0;

        if (!std::isfinite(upper * upper))
        {throw std::runtime_error("Unable to find a Mach upper bound.");}
    }

    // narrow the bounds
    for (int iteration = 0; iteration < 100; ++iteration)
    {
        const double midpoint = lower + (upper - lower) / 2.0;

        if (midpoint == lower || midpoint == upper)
        {break;}

        if (prandtlMeyer(midpoint, gamma) < nu)
        {lower = midpoint;}
        else
        {upper = midpoint;}
    }

    return lower + (upper - lower) / 2.0;
}

int main()
{
    try
    {
        // flow inputs
        const double gamma = 1.4;
        const double inletMach = 1.775;
        const double turnDegrees = 16.0;

        // static inlet conditions in pa and kelvin
        const double inletPressure = 101325.0;
        const double inletTemperature = 273.15;

        // geometry inputs
        const double cornerX = 0.0;
        const double cornerY = 0.0;
        const double rayLength = 10.0;

        // one segment gives two endpoints per ray
        const int intervals = 100000;
        const int segmentsPerRay = 1;

        // check inputs
        if (intervals <= 0 || segmentsPerRay <= 0)
        {throw std::invalid_argument("Intervals and segments must be positive.");}

        if (!std::isfinite(turnDegrees) || turnDegrees < 0.0 ||
            !std::isfinite(rayLength) || rayLength <= 0.0 ||
            !std::isfinite(cornerX) || !std::isfinite(cornerY))
        {throw std::invalid_argument("Invalid geometry inputs.");}

        if (!std::isfinite(inletPressure) || inletPressure <= 0.0 ||
            !std::isfinite(inletTemperature) || inletTemperature <= 0.0)
        {throw std::invalid_argument("Pressure and temperature must be positive.");}

        const double totalTurn = radians(turnDegrees);
        const double inletNu = prandtlMeyer(inletMach, gamma);

        // check the final expansion before starting threads
        const double outletMach =
            machFromNu(inletNu + totalTurn, gamma);

        const double inletFactor =
            1.0 + 0.5 * (gamma - 1.0) * inletMach * inletMach;

        std::vector<FanState> states(static_cast<std::size_t>(intervals) + 1);

        const double calculationStart = omp_get_wtime();

        // calculate each ray independently
        #pragma omp parallel for schedule(static)
        for (int i = 0; i <= intervals; ++i)
        {
            const double fraction =
                static_cast<double>(i) / intervals;

            const double delta = totalTurn * fraction;

            FanState state{};

            state.theta = -delta;
            state.nu = inletNu + delta;
            state.Mach = machFromNu(state.nu, gamma);
            state.mu = std::asin(1.0 / state.Mach);
            state.alpha = state.theta + state.mu;

            const double localFactor =
                1.0 + 0.5 * (gamma - 1.0) *
                state.Mach * state.Mach;

            const double temperatureRatio =
                inletFactor / localFactor;

            state.temperature =
                inletTemperature * temperatureRatio;

            state.pressure =
                inletPressure *
                std::pow(
                    temperatureRatio,
                    gamma / (gamma - 1.0)
                );

            states[i] = state;
        }

        const double calculationTime =
            omp_get_wtime() - calculationStart;

        const double exportStart = omp_get_wtime();

        // save in the current working directory
        const auto outputPath =
            std::filesystem::absolute("characteristics.csv");

        std::ofstream csv(outputPath);

        if (!csv)
        {throw std::runtime_error("Cannot open CSV output.");}

        csv << std::setprecision(15);

        csv << "line_id,point_id,x,y,Mach,"
               "theta_deg,nu_deg,mu_deg,alpha_deg,"
               "pressure_Pa,temperature_K\n";

        // write the results in order
        for (int i = 0; i <= intervals; ++i)
        {
            const FanState& state = states[i];

            const double directionX = std::cos(state.alpha);
            const double directionY = std::sin(state.alpha);

            const double thetaDegrees = degrees(state.theta);
            const double nuDegrees = degrees(state.nu);
            const double muDegrees = degrees(state.mu);
            const double alphaDegrees = degrees(state.alpha);

            for (int j = 0; j <= segmentsPerRay; ++j)
            {
                const double s =
                    rayLength *
                    (static_cast<double>(j) / segmentsPerRay);

                const double x = cornerX + s * directionX;
                const double y = cornerY + s * directionY;

                csv << i << ',' << j << ','
                    << x << ',' << y << ','
                    << state.Mach << ','
                    << thetaDegrees << ','
                    << nuDegrees << ','
                    << muDegrees << ','
                    << alphaDegrees << ','
                    << state.pressure << ','
                    << state.temperature << '\n';
            }
        }

        csv.close();

        if (!csv)
        {throw std::runtime_error("CSV writing failed.");}

        const double exportTime = omp_get_wtime() - exportStart;

        // print a short summary
        std::cout << std::fixed << std::setprecision(4);

        std::cout << "Characteristics: " << states.size() << '\n';
        std::cout << "Outlet Mach: " << outletMach << '\n';
        std::cout << "Outlet pressure: "
                  << states.back().pressure / 1000.0 << " kPa\n";
        std::cout << "Outlet temperature: "
                  << states.back().temperature << " K\n";
        std::cout << "Calculation time: "
                  << calculationTime << " s\n";
        std::cout << "CSV export time: "
                  << exportTime << " s\n";
        std::cout << "Saved: " << outputPath.string() << '\n';
    }
    catch (const std::exception& error)
    {std::cerr << "Error: " << error.what() << '\n';
    return 1;}

    return 0;
}